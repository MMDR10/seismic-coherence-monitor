# 🌏 Seismic Coherence Monitor

每日自動監測台灣近場站 spectral coherence — D₁ 框架 cross-domain 應用（火山×地震耦合測試延伸）。

## 背景

2026-08-01 火山×地震耦合測試發現：**spectral coherence (0.5–5 Hz) 係地震域可靠 discriminator**（火山 dH_curl 方法對地震無效，但 coherence 捕捉到波前相干）：

| 事件 | bg | cos | Δ | p |
|------|-----|-----|-----|-----|
| 2022 關山 M6.4 | 0.065 | 0.139 | +0.074 | **4.5e-17** |
| 2025 嘉義 M6.0 | 0.065 | 0.142 | +0.077 | **1.9e-19** |
| 2018 花蓮 M6.4 | 0.080 | 0.109 | +0.028 | **1.6e-6** |

- 地震波前到達時全網 coherence 平均抬升 ~2×，震後 1–6h 持續（p<1e-48）
- 控制日 mean=0.068、p99=0.205 — 有零星高值，判據要「分佈整體抬升」而非絕對門檻

## 監測內容（每日 05:30 UTC 自動跑）

1. **FDSN 拉取**：TW 網絡近場站 NACB/YHNB/SSLB 最近 24h 波形（20 Hz）
2. **計算**：0.5–5 Hz 帶通 → 120s/60s sliding window 全網平均 spectral coherence
3. **USGS 交叉驗證**：當日台灣 M≥4.5 地震目錄
4. **判據**：對比控制日 baseline p99=0.205（>p99 比例 + 最長連續窗口）
5. **記錄**：append 至 `seismic_history.json`（長期數據收集）

## v2 升級（2026-09-01）

回應 8/25 M5.5（Hengchun NE）個案：v1「分佈抬升」判據漏報瞬時 spike（max_coh=0.717 @ 地震前 1 分鐘，但 frac>p99 僅 0.5%）。v2 新增：

- **spike 偵測（L1）**：`max_coh` 時刻 + robust z-score（vs 全日 MAD）+ max/median ratio；與 USGS 事件 ±30 分鐘對齊檢查 → `spike_alert`
- **時間序列存檔（L1）**：`coh_series`（全日 120s/60s 值）寫入每日 result，之後唔使重跑先拎到峰值時刻
- **雙頻帶對照（L2）**：新增 0.1–0.5 Hz band（火山研究實錘 lock 帶），含第二頻帶 coherence + Hilbert 相位鎖定 `c=cos(Δφ)` kappa
- **判據並行**：`alert`（v1 分佈抬升型）與 `spike_alert`（v2 瞬時+事件對齊型）獨立輸出，歷史 schema 向後相容

驗證（本機重算）：8/25 M5.5 → `spike_alert=True`（z=89.7，峰值 06:59 UTC = 地震前 1 分鐘）；8/15 M4.9 → `spike_alert=True`（峰值 11:29 UTC = 地震前 1 分鐘）；8/01 無地震 → `spike_alert=False`（無誤報）。兩次事件峰值均出現於地震前 1 分鐘，pattern 一致。

## 使用

```bash
# 指定日期（前一日 24h）
python seismic_monitor.py --date 20260801

# 最近 24h（預設）
python seismic_monitor.py

# 自訂窗口 / 站
python seismic_monitor.py --hours 6 --stations NACB,YHNB,SSLB
```

## 輸出

- `seismic_result_YYYYMMDD.json` — 當日測量
- `seismic_history.json` — 累積時間序列

## 警報邏輯（v2 雙警報）

**`alert=True`（v1 分佈抬升型，對應火山 tremor 式持續抬升）**：
>baseline p99 比例 ≥5% **且** 最長連續 ≥8 窗口（=8+ 分鐘持續抬升）。
有 M≥5 地震時高 coherence 係預期行為；無地震時高 coherence 先係異常（待查）。

**`spike_alert=True`（v2 瞬時+對齊型，對應地震 spike）**：
全日 max_coh 為 spike（robust z ≥ 6 或 max_coh ≥ 0.55）**且** 與 USGS M≥4.5 事件 ±30 分鐘對齊。
v1 對瞬時 spike 免疫（8/25 M5.5 個案），v2 補返呢個盲區。

## 作者

tygtDc, Deep Research（D₁ 框架跨域研究）
