#!/usr/bin/env python3
"""
Seismic Coherence Monitor v3 — 每日自動監測台灣近場站 spectral coherence + spike 偵測
======================================================================================
v2 升級（2026-09-01, DR 補足 — 回應 8/26 M5.5 個案盲區）：
  L1. spike 偵測：max_coh 時刻 + max/median ratio + robust z-score；USGS event 時間對齊
  L1. 存完整 coherence 時間序列（coh_series）→ 之後唔使重跑先拎到峰值時刻
  L2. 第二頻帶 0.1-0.5Hz coherence + 相位鎖定 kappa（c=cos(Δφ)，3 站 pair 平均）
      —— 火山研究實錘 broad 0.5-5Hz 係「反轉帶」，真正鎖定喺 0.1-0.5Hz；做雙頻帶對照
  v1 判據保留（>p99 比例 ≥5% / 連續 ≥8 窗，分佈抬升型），新增 spike_alert（瞬時型）

v3 升級（2026-11-13, DR — 10 無震日重校 + rolling 自校準）：
  L1. 常數重校：CONTROL_P99 0.205→0.145 / CONTROL_MEAN 0.068→0.066 / ABS_SPIKE 0.55→0.31
      （舊 2022-09-14 單控制日 baseline 已 stale：ABS_SPIKE 跌入無震日背景 0.30–0.58；
      10 無震日 23,510 窗 pooled 校準，FP 驗證 0/10，見 REPORT3）
  L2. rolling 30-day 自校準：由 seismic_history.json 每日摘要 median-of-daily-p99 自動追蹤
      背景漂移；污染防護 = 排除事件日（spike_alert / n_eq_ge50>0）+ clamp [0.5×,2.5×]；
      乾淨日 <7 → fallback 固定常數。記錄 baseline_source 於 result
  L3. v1 連續窗門檻 8→10（--run-threshold）：10 無震日 max run=5，margin 1.6×→2×
  L4. Z_SPIKE 標記為裝飾性：無震日 in-day z 39–83 恒過線（重尾），v2 真 gate 係事件對齊

基於 D₁ 框架 cross-domain 發現（2026-08-01 耦合測試）：spectral coherence 係地震域
reliable discriminator；控制日 (2022-09-14) mean=0.068, p99=0.205, max=0.392。

每日任務：
  1. 拉取 TW 網絡近場站（NACB/YHNB/SSLB）最近 24h 連續波形
  2. 0.5-5Hz 帶通 → 120s/60s sliding window spectral coherence（全網平均）
  3. 0.1-0.5Hz 帶通 → 同窗 coherence + Hilbert 相位鎖定 kappa（雙頻帶對照）
  4. USGS catalog 交叉驗證：當日有冇 M≥4.5 台灣地震；max_coh 時刻 ±30min 對齊檢查
  5. 判據：v1 分佈抬升（>p99 比例/連續窗口）+ v2 spike（robust z / 絕對閾值 × 事件對齊）
  6. Append 到 seismic_history.json + 生成當日報告（含 coh_series 時間序列）

用法：
  python seismic_monitor.py [--date YYYYMMDD] [--hours N] [--baseline-p99 0.205]
                            [--stations NACB,YHNB,SSLB] [--outdir .]

輸出：
  seismic_result_YYYYMMDD.json — 當日測量（含 spike / event_alignment / coh_series）
  seismic_history.json         — 累積時間序列（含 spike 摘要）
"""
import argparse, json, os, sys, time, warnings
import numpy as np
from scipy.signal import butter, sosfiltfilt, coherence, hilbert
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import urllib.request

warnings.filterwarnings('ignore')

NETWORK = 'TW'
DEFAULT_STATIONS = ['NACB', 'YHNB', 'SSLB']
BAND_LO, BAND_HI = 0.5, 5.0          # v1 主頻帶
BAND2_LO, BAND2_HI = 0.1, 0.5        # v2 第二頻帶（火山研究鎖定帶）
WIN_S, STEP_S = 120, 60
NPERSEG = 256
# ---- v3 重校常數（2026-11-13；10 無震日 pooled 23,510 窗：TW 10-02/08/14/19/25/31, 11-05/07/10/12）----
CONTROL_P99 = 0.145   # pooled p99（舊 0.205=2022-09-14 單日，偏高損靈敏度；見 REPORT3）
CONTROL_MEAN = 0.066  # pooled mean（舊 0.068，實質唔變）
BAND2_P99_REF = 0.133 # 0.1-0.5Hz coherence 參考線（pooled p99）
KAPPA_P99_REF = 0.183 # 0.1-0.5Hz kappa 參考線（pooled p99）
# v2 spike 閾值（v3 重校）
Z_SPIKE = 6.0          # ⚠️ 裝飾性：無震日 in-day z 39–83 恒過線（重尾），實質永不 binding
ABS_SPIKE = 0.31       # optB = pooled p99.9（舊 0.55 跌入背景 0.30–0.58；0.31 保留 M3.x 同震偵測）
ALIGN_MIN = 30         # max_coh 時刻與 M≥4.5 事件 ±30 分鐘
# ---- v3 rolling 自校準 ----
RUN_THRESHOLD = 10     # v1 連續窗門檻（v2=8；10 無震日 max run=5 → margin 2×）
ROLLING_DAYS = 30      # rolling 窗口（日數）
ROLLING_MIN_DAYS = 7   # 乾淨日 < 呢個數 → fallback 固定常數
ROLLING_CLAMP = (0.5, 2.5)  # rolling p99 相對 CONTROL_P99 嘅 clamp（防 history 污染跑飛）
USGS_URL = ('https://earthquake.usgs.gov/fdsnws/event/1/query'
            '?format=geojson&starttime={}&endtime={}&minmagnitude=4.5'
            '&minlatitude=21.5&maxlatitude=25.5&minlongitude=119.5&maxlongitude=122.5')


def bandpass(x, fs, lo=BAND_LO, hi=BAND_HI):
    sos = butter(4, [lo, hi], btype='band', fs=fs, output='sos')
    return sosfiltfilt(sos, x)


def spec_coh_pair(x, y, fs, lo=BAND_LO, hi=BAND_HI, nperseg=NPERSEG):
    try:
        f, C = coherence(x, y, fs=fs, nperseg=nperseg, noverlap=nperseg // 2)
        mask = (f >= lo) & (f <= hi)
        if np.sum(mask) > 3:
            return float(np.mean(C[mask]))
    except Exception:
        pass
    return 0.0


def phase_lock_kappa_pair(x, y, fs, lo=BAND2_LO, hi=BAND2_HI):
    """Hilbert 即時相位差集中度 c=cos(Δφ)（0.1-0.5Hz 鎖定帶）"""
    try:
        xb = bandpass(x, fs, lo, hi)
        yb = bandpass(y, fs, lo, hi)
        hx = hilbert(xb)
        hy = hilbert(yb)
        dphi = np.angle(hx * np.conj(hy))
        return float(np.mean(np.cos(dphi)))
    except Exception:
        return 0.0


def fetch_waveforms(client, stations, start, end, attempts=3):
    """逐站拉取，pair-wise 容錯（唔強制共同窗口）"""
    data_list, srs = [], []
    for st in stations:
        got = False
        for attempt in range(attempts):
            try:
                tr = client.get_waveforms(
                    NETWORK, st, '*', 'HH?,BH?,EH?', start, end,
                    attach_response=False)[0]
                data_list.append(tr.data.astype(np.float64))
                srs.append(tr.stats.sampling_rate)
                got = True
                break
            except Exception as e:
                time.sleep(2 * (attempt + 1))
        if not got:
            print(f'  ⚠️ {st}: 拉取失敗 (skip)')
    return data_list, srs


def daily_series(data_list, srs, lo, hi, metric='coh'):
    """全日 120s/60s sliding window 全網平均 spectral coherence / kappa"""
    Nmin = min(len(x) for x in data_list)
    fs = min(srs)
    win = int(WIN_S * fs)
    step = int(STEP_S * fs)
    ts, vals = [], []
    i = 0
    while i + win <= Nmin:
        cs = []
        for a in range(len(data_list)):
            for b in range(a + 1, len(data_list)):
                if metric == 'coh':
                    cs.append(spec_coh_pair(data_list[a][i:i + win],
                                            data_list[b][i:i + win], fs, lo, hi))
                else:  # kappa
                    cs.append(phase_lock_kappa_pair(data_list[a][i:i + win],
                                                    data_list[b][i:i + win], fs, lo, hi))
        vals.append(float(np.mean(cs)))
        ts.append(i / fs)  # seconds since window start
        i += step
    return np.array(ts), np.array(vals)


def usgs_events(start, end):
    """當日台灣 M≥4.5 地震目錄（交叉驗證用）"""
    url = USGS_URL.format(start.strftime('%Y-%m-%dT%H:%M:%S'),
                          end.strftime('%Y-%m-%dT%H:%M:%S'))
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.load(r)
        evs = []
        for f in data.get('features', []):
            p = f['properties']
            g = f['geometry']['coordinates']
            evs.append({'time': p.get('time'), 'mag': p.get('mag'),
                        'lon': g[0], 'lat': g[1], 'place': p.get('place', '')})
        evs.sort(key=lambda e: e['time'])
        return evs
    except Exception as e:
        print(f'  ⚠️ USGS query 失敗: {str(e)[:120]}')
        return None


def longest_run(series):
    best = cur = 0
    for v in series:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


def robust_z(x, med, mad):
    """robust z-score: (x - median) / (1.4826 * MAD)"""
    if mad <= 0 or not np.isfinite(mad):
        return 0.0
    return float((x - med) / (1.4826 * mad))


def detect_spike(ts, cohs, start, evs):
    """v2 spike 偵測：max 時刻 + robust z + event 對齊（±30min）"""
    out = {'max_coh': float(np.max(cohs)),
           'max_coh_time_utc': str(start + ts[int(np.argmax(cohs))]),
           'median_coh': float(np.median(cohs)),
           'p99_coh': float(np.percentile(cohs, 99))}
    med, mad = float(np.median(cohs)), float(np.median(np.abs(cohs - np.median(cohs))))
    out['mad'] = mad
    out['max_median_ratio'] = float(np.max(cohs) / med) if med > 0 else 0.0
    out['max_robust_z'] = robust_z(out['max_coh'], med, mad)
    # top-5 spikes
    idx = np.argsort(cohs)[-5:][::-1]
    out['top_spikes'] = [{'t_utc': str(start + ts[j]),
                          'coh': float(cohs[j])} for j in sorted(idx)]
    # event 對齊
    peak_t = start + ts[int(np.argmax(cohs))]
    aligned = None
    if evs:
        for e in evs:
            et = UTCDateTime(e['time'] / 1000)
            dmin = abs((peak_t - et) / 60.0)
            if dmin <= ALIGN_MIN:
                aligned = {'event_time_utc': str(et), 'mag': e['mag'],
                           'place': e['place'], 'delta_min': round(float(dmin), 1)}
                break
    out['event_aligned'] = aligned is not None
    out['aligned_event'] = aligned
    # spike 判定：robust z 高 或 max 絕對高，且對齊事件
    spike_strong = (out['max_robust_z'] >= Z_SPIKE or out['max_coh'] >= ABS_SPIKE)
    out['spike_alert'] = bool(spike_strong and out['event_aligned'])
    out['spike_strong_no_event'] = bool(spike_strong and not out['event_aligned'])
    return out


def rolling_baseline(hist, days=ROLLING_DAYS, min_days=ROLLING_MIN_DAYS):
    """v3 rolling 自校準：由 history records 嘅每日摘要（p99_coh/mean_coh）估 rolling baseline。
    estimator = median-of-daily-p99（對單日污染 robust；10 無震日驗證 ≈ pooled p99）。
    污染防護：① 排除事件日（spike_alert=True 或 n_eq_ge50>0）② clamp [0.5×,2.5×]×CONTROL_P99。
    回傳 dict；唔夠乾淨日回傳 None（caller fallback 固定常數）。
    註：history schema 只存每日摘要（唔存全序列），所以唔可以 pooled windows；
    band B / kappa 暫用固定參考線（BAND2_P99_REF / KAPPA_P99_REF）。
    """
    recs = hist.get('records', []) if isinstance(hist, dict) else []
    if not recs:
        return None
    recs = sorted(recs, key=lambda r: r.get('date', ''))[-days:]
    clean = [r for r in recs
             if not r.get('spike_alert') and not (r.get('n_eq_ge50') or 0) > 0
             and r.get('p99_coh') is not None]
    if len(clean) < min_days:
        return None
    p99s = np.array([r['p99_coh'] for r in clean], dtype=float)
    means = np.array([r.get('mean_coh') for r in clean if r.get('mean_coh') is not None], dtype=float)
    med = float(np.median(p99s))
    lo, hi = ROLLING_CLAMP[0] * CONTROL_P99, ROLLING_CLAMP[1] * CONTROL_P99
    p99 = float(np.clip(med, lo, hi))
    return {
        'p99': round(p99, 4),
        'mean': round(float(np.median(means)), 4) if means.size else None,
        'n_days_used': len(clean),
        'n_days_excluded': len(recs) - len(clean),
        'daily_p99_range': [round(float(p99s.min()), 4), round(float(p99s.max()), 4)],
        'clamped': bool(med < lo or med > hi),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--date', default=None, help='YYYYMMDD (default: yesterday UTC)')
    ap.add_argument('--hours', type=int, default=24)
    ap.add_argument('--stations', default=','.join(DEFAULT_STATIONS))
    ap.add_argument('--baseline-p99', type=float, default=CONTROL_P99,
                    help='固定 fallback p99（rolling 數據唔夠時用）')
    ap.add_argument('--run-threshold', type=int, default=RUN_THRESHOLD,
                    help='v1 連續超標窗口門檻（v3: 10）')
    ap.add_argument('--outdir', default='.')
    args = ap.parse_args()

    stations = [s.strip() for s in args.stations.split(',') if s.strip()]
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # v3 rolling 自校準（先讀 history 決定 baseline，再做當日分析）
    hpath = os.path.join(outdir, 'seismic_history.json')
    hist_existing = None
    if os.path.exists(hpath):
        with open(hpath) as f:
            hist_existing = json.load(f)
    roll = rolling_baseline(hist_existing) if hist_existing else None
    if roll:
        baseline_p99 = roll['p99']
        baseline_source = f'rolling_{ROLLING_DAYS}d(n={roll["n_days_used"]})'
    else:
        baseline_p99 = args.baseline_p99
        baseline_source = f'fixed({args.baseline_p99})'

    if args.date:
        end = UTCDateTime(args.date + 'T00:00:00')
        start = end - args.hours * 3600
        day_tag = args.date
    else:
        end = UTCDateTime() - 3600  # 1h latency for data availability
        start = end - args.hours * 3600
        day_tag = end.strftime('%Y%m%d')

    print(f'Seismic Coherence Monitor v3 — {day_tag}')
    print(f'  window: {start} → {end} ({args.hours}h)')
    print(f'  stations: {stations}  baseline p99={baseline_p99:.4f} [{baseline_source}]')

    # 1) USGS 交叉驗證（先查有冇地震，再決定睇法）
    evs = usgs_events(start, end)
    big = []
    if evs is not None:
        big = [e for e in evs if (e['mag'] or 0) >= 5.0]
        print(f'  USGS: {len(evs)} 個 M≥4.5 事件, {len(big)} 個 M≥5.0')
        for e in evs[:5]:
            t = UTCDateTime(e['time'] / 1000)
            print(f'    {t.strftime("%m-%d %H:%M")} M{e["mag"]} @ '
                  f'({e["lat"]:.2f},{e["lon"]:.2f}) {e["place"][:50]}')

    # 2) FDSN 拉波形
    client = Client('EARTHSCOPE')
    data_list, srs = fetch_waveforms(client, stations, start, end)
    if len(data_list) < 2:
        print('❌ 少於 2 站數據，無法計算 coherence')
        sys.exit(1)
    print(f'  拉取成功: {len(data_list)} 站, 採樣率 {srs}')

    # 3) 全日 coherence series（0.5-5Hz 主帶）
    ts, cohs = daily_series(data_list, srs, BAND_LO, BAND_HI, metric='coh')
    print(f'  全日窗口: {len(cohs)} 個 (120s/60s step)')
    print(f'  [0.5-5Hz] mean={np.mean(cohs):.4f}  max={np.max(cohs):.4f}  '
          f'p95={np.percentile(cohs, 95):.4f}  p99={np.percentile(cohs, 99):.4f}')

    # 3b) 第二頻帶 0.1-0.5Hz coherence + kappa（火山研究鎖定帶對照）
    ts2, cohs2 = daily_series(data_list, srs, BAND2_LO, BAND2_HI, metric='coh')
    ts_k, kap = daily_series(data_list, srs, BAND2_LO, BAND2_HI, metric='kappa')
    print(f'  [0.1-0.5Hz] mean={np.mean(cohs2):.4f}  max={np.max(cohs2):.4f}  '
          f'p99={np.percentile(cohs2, 99):.4f}  kappa mean={np.mean(kap):.4f} '
          f'max={np.max(kap):.4f}')

    # 4) v1 判據 vs baseline（分佈抬升型；v3: rolling baseline + run 門檻 10）
    frac_above = float(np.mean(cohs > baseline_p99))
    run = longest_run(cohs > baseline_p99)
    mean_ratio = float(np.mean(cohs) / CONTROL_MEAN) if CONTROL_MEAN else 0.0
    alert = (frac_above >= 0.05 or run >= args.run_threshold)
    if evs is not None and big:
        alert = alert and True  # 有地震時高 coherence 係預期行為
    else:
        alert = alert and False  # 無地震時高 coherence 先係異常（待查）

    # 4b) v2 spike 偵測（瞬時型）
    spike = detect_spike(ts, cohs, start, evs)

    print(f'  >baseline p99 ({baseline_p99:.4f} @ {baseline_source}): {100*frac_above:.1f}%  '
          f'最長連續 {run} 窗口 (門檻 {args.run_threshold})')
    print(f'  mean / 控制日 mean ratio: {mean_ratio:.2f}x')
    print(f'  [v1] 🚨 ALERT(分佈抬升): {alert}')
    print(f'  [v2] spike: max={spike["max_coh"]:.4f} @ {spike["max_coh_time_utc"][11:16]}UTC '
          f'z={spike["max_robust_z"]:.1f} ratio={spike["max_median_ratio"]:.1f}x '
          f'aligned={spike["event_aligned"]}')
    print(f'  [v2] 🚨 SPIKE_ALERT(瞬時+對齊): {spike["spike_alert"]}')

    # 5) 記錄
    result = {
        'date': day_tag,
        'window_start': str(start), 'window_end': str(end),
        'stations': stations,
        'n_windows': len(cohs),
        'mean_coh': float(np.mean(cohs)),
        'max_coh': float(np.max(cohs)),
        'p95_coh': float(np.percentile(cohs, 95)),
        'p99_coh': float(np.percentile(cohs, 99)),
        'baseline_p99': baseline_p99,
        'baseline_source': baseline_source,
        'rolling': roll,
        'run_threshold': args.run_threshold,
        'frac_above_p99': frac_above,
        'longest_run_above_p99': int(run),
        'mean_ratio_vs_control': mean_ratio,
        'n_eq_ge45': len(evs) if evs is not None else None,
        'n_eq_ge50': len(big) if evs is not None else None,
        'earthquakes': evs if evs else [],
        'alert': bool(alert),
        # v2 新增
        'spike': spike,
        'spike_alert': spike['spike_alert'],
        'band_01_05': {
            'mean_coh': float(np.mean(cohs2)),
            'max_coh': float(np.max(cohs2)),
            'p99_coh': float(np.percentile(cohs2, 99)),
            'max_kappa': float(np.max(kap)),
            'mean_kappa': float(np.mean(kap)),
            'max_coh_time_utc': str(start + ts2[int(np.argmax(cohs2))]),
        },
        'coh_series': {
            't_sec': [int(x) for x in ts.tolist()],
            'coh': [round(float(x), 4) for x in cohs.tolist()],
        },
        'updated_utc': UTCDateTime().strftime('%Y-%m-%dT%H:%M:%S'),
    }
    rpath = os.path.join(outdir, f'seismic_result_{day_tag}.json')
    with open(rpath, 'w') as f:
        json.dump(result, f, indent=2)

    # 6) Append history（數據收集，保留 v1 schema + v2 spike 摘要）
    hpath = os.path.join(outdir, 'seismic_history.json')
    if os.path.exists(hpath):
        with open(hpath) as f:
            hist = json.load(f)
    else:
        hist = {'control_baseline': {'mean': CONTROL_MEAN, 'p99': CONTROL_P99,
                                     'desc': '2022-09-14 控制日 TW NACB/YHNB/SSLB'},
                'records': []}
    hist['records'] = [r for r in hist['records'] if r.get('date') != day_tag]
    hist['records'].append({k: result[k] for k in
                            ['date', 'mean_coh', 'max_coh', 'p99_coh',
                             'frac_above_p99', 'longest_run_above_p99',
                             'n_eq_ge50', 'alert']} | {
        'spike_alert': spike['spike_alert'],
        'spike_max_coh': spike['max_coh'],
        'spike_max_z': spike['max_robust_z'],
        'spike_time_utc': spike['max_coh_time_utc'],
        'spike_aligned': spike['event_aligned'],
        'band01_05_max': float(np.max(cohs2)),
        'band01_05_kappa_max': float(np.max(kap)),
        'band01_05_p99': float(np.percentile(cohs2, 99)),
        'band01_05_kappa_p99': float(np.percentile(kap, 99)),
    })
    hist['records'].sort(key=lambda r: r['date'])
    with open(hpath, 'w') as f:
        json.dump(hist, f, indent=2)
    print(f'✅ 完成: {rpath}')
    print(f'   歷史累積: {len(hist["records"])} 日')


if __name__ == '__main__':
    main()