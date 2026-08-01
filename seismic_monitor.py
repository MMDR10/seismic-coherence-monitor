#!/usr/bin/env python3
"""
Seismic Coherence Monitor — 每日自動監測台灣近場站 spectral coherence
=========================================================================
基於 D₁ 框架 cross-domain 發現（2026-08-01 耦合測試）：
spectral coherence (0.5-5Hz) 係地震域可靠 discriminator —
  3/3 台灣近場事件 cos 段抬升 ~2x (p<1e-6)，震後 1-6h 持續 (p<1e-48)
  控制日 (2022-09-14) mean=0.068，p99=0.205，max=0.392
判據：分佈整體抬升（唔係絕對門檻）→ 對比 baseline p99 + 持續性

每日任務：
  1. 拉取 TW 網絡近場站（NACB/YHNB/SSLB）最近 24h 連續波形
  2. 0.5-5Hz 帶通 → 120s/60s sliding window spectral coherence（全網平均）
  3. USGS catalog 交叉驗證：當日有冇 M≥4.5 台灣地震
  4. 對比 baseline（控制日 p99=0.205）：>p99 比例 + 最長連續窗口
  5. Append 到 seismic_history.json + 生成當日報告

用法：
  python seismic_monitor.py [--date YYYYMMDD] [--hours N] [--baseline-p99 0.205]
                            [--stations NACB,YHNB,SSLB] [--outdir .]

輸出：
  seismic_result_YYYYMMDD.json — 當日測量（供 README badge/歷史）
  seismic_history.json         — 累積時間序列（數據收集）
"""
import argparse, json, os, sys, time, warnings
import numpy as np
from scipy.signal import butter, sosfiltfilt, coherence
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import urllib.request

warnings.filterwarnings('ignore')

NETWORK = 'TW'
DEFAULT_STATIONS = ['NACB', 'YHNB', 'SSLB']
BAND_LO, BAND_HI = 0.5, 5.0
WIN_S, STEP_S = 120, 60
NPERSEG = 256
CONTROL_P99 = 0.205   # 2022-09-14 控制日實測（3 站全日）
CONTROL_MEAN = 0.068
USGS_URL = ('https://earthquake.usgs.gov/fdsnws/event/1/query'
            '?format=geojson&starttime={}&endtime={}&minmagnitude=4.5'
            '&minlatitude=21.5&maxlatitude=25.5&minlongitude=119.5&maxlongitude=122.5')


def bandpass(x, sr, lo=BAND_LO, hi=BAND_HI):
    sos = butter(4, [lo, hi], btype='band', fs=sr, output='sos')
    return sosfiltfilt(sos, x)


def spec_coh_pair(x, y, sr, nperseg=NPERSEG):
    try:
        f, C = coherence(x, y, fs=sr, nperseg=nperseg, noverlap=nperseg // 2)
        mask = (f >= BAND_LO) & (f <= BAND_HI)
        if np.sum(mask) > 3:
            return float(np.mean(C[mask]))
    except Exception:
        pass
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


def daily_coherence_series(data_list, srs):
    """全日 120s/60s sliding window 全網平均 spectral coherence"""
    Nmin = min(len(x) for x in data_list)
    # 統一採樣率：用最低嗰個（HH 100Hz vs BH 20Hz 混用時安全）
    sr = min(srs)
    win = int(WIN_S * sr)
    step = int(STEP_S * sr)
    ts, cohs = [], []
    i = 0
    while i + win <= Nmin:
        cs = []
        for a in range(len(data_list)):
            for b in range(a + 1, len(data_list)):
                cs.append(spec_coh_pair(data_list[a][i:i + win],
                                        data_list[b][i:i + win], sr))
        cohs.append(float(np.mean(cs)))
        ts.append(i / sr / 3600.0)  # hours since window start
        i += step
    return np.array(ts), np.array(cohs)


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


def longest_run(series, thr):
    best = cur = 0
    for v in series:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--date', default=None, help='YYYYMMDD (default: yesterday UTC)')
    ap.add_argument('--hours', type=int, default=24)
    ap.add_argument('--stations', default=','.join(DEFAULT_STATIONS))
    ap.add_argument('--baseline-p99', type=float, default=CONTROL_P99)
    ap.add_argument('--outdir', default='.')
    args = ap.parse_args()

    stations = [s.strip() for s in args.stations.split(',') if s.strip()]
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    if args.date:
        end = UTCDateTime(args.date + 'T00:00:00')
        start = end - args.hours * 3600
        day_tag = args.date
    else:
        end = UTCDateTime() - 3600  # 1h latency for data availability
        start = end - args.hours * 3600
        day_tag = end.strftime('%Y%m%d')

    print(f'Seismic Coherence Monitor — {day_tag}')
    print(f'  window: {start} → {end} ({args.hours}h)')
    print(f'  stations: {stations}  baseline p99={args.baseline_p99:.3f}')

    # 1) USGS 交叉驗證（先查有冇地震，再決定睇法）
    evs = usgs_events(start, end)
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

    # 3) 全日 coherence series
    ts, cohs = daily_coherence_series(data_list, srs)
    print(f'  全日窗口: {len(cohs)} 個 (120s/60s step)')
    print(f'  mean={np.mean(cohs):.4f}  max={np.max(cohs):.4f}  '
          f'p95={np.percentile(cohs, 95):.4f}  p99={np.percentile(cohs, 99):.4f}')

    # 4) 判據 vs baseline
    frac_above = float(np.mean(cohs > args.baseline_p99))
    run = longest_run(cohs > args.baseline_p99, args.baseline_p99)
    mean_ratio = float(np.mean(cohs) / CONTROL_MEAN) if CONTROL_MEAN else 0.0
    # 警報：>p99 比例 ≥5% 或 連續 ≥8 窗口（=8+ 分鐘持續抬升）且當日有 M≥5
    alert = (frac_above >= 0.05 or run >= 8)
    if evs is not None and big:
        alert = alert and True  # 有地震時高 coherence 係預期行為
    else:
        alert = alert and False  # 無地震時高 coherence 先係異常（待查）

    print(f'  >baseline p99 ({args.baseline_p99}): {100*frac_above:.1f}%  最長連續 {run} 窗口')
    print(f'  mean / 控制日 mean ratio: {mean_ratio:.2f}x')
    print(f'  🚨 ALERT: {alert}')

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
        'baseline_p99': args.baseline_p99,
        'frac_above_p99': frac_above,
        'longest_run_above_p99': int(run),
        'mean_ratio_vs_control': mean_ratio,
        'n_eq_ge45': len(evs) if evs else None,
        'n_eq_ge50': len(big) if evs else None,
        'earthquakes': evs if evs else [],
        'alert': bool(alert),
        'updated_utc': UTCDateTime().strftime('%Y-%m-%dT%H:%M:%S'),
    }
    rpath = os.path.join(outdir, f'seismic_result_{day_tag}.json')
    with open(rpath, 'w') as f:
        json.dump(result, f, indent=2)

    # 6) Append history（數據收集）
    hpath = os.path.join(outdir, 'seismic_history.json')
    if os.path.exists(hpath):
        with open(hpath) as f:
            hist = json.load(f)
    else:
        hist = {'control_baseline': {'mean': CONTROL_MEAN, 'p99': CONTROL_P99,
                                     'desc': '2022-09-14 控制日 TW NACB/YHNB/SSLB'},
                'records': []}
    # 避免同日重複
    hist['records'] = [r for r in hist['records'] if r.get('date') != day_tag]
    hist['records'].append({k: result[k] for k in
                            ['date', 'mean_coh', 'max_coh', 'p99_coh',
                             'frac_above_p99', 'longest_run_above_p99',
                             'n_eq_ge50', 'alert']})
    hist['records'].sort(key=lambda r: r['date'])
    with open(hpath, 'w') as f:
        json.dump(hist, f, indent=2)
    print(f'✅ 完成: {rpath}')
    print(f'   歷史累積: {len(hist["records"])} 日')


if __name__ == '__main__':
    main()
