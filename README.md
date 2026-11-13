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
4. **判據**：對比 baseline p99（v3：rolling 30-day 自校準，fallback 0.145；>p99 比例 + 最長連續窗口）
5. **記錄**：append 至 `seismic_history.json`（長期數據收集）

## v2 升級（2026-09-01）

回應 8/25 M5.5（Hengchun NE）個案：v1「分佈抬升」判據漏報瞬時 spike（max_coh=0.717 @ 地震前 1 分鐘，但 frac>p99 僅 0.5%）。v2 新增：

- **spike 偵測（L1）**：`max_coh` 時刻 + robust z-score（vs 全日 MAD）+ max/median ratio；與 USGS 事件 ±30 分鐘對齊檢查 → `spike_alert`
- **時間序列存檔（L1）**：`coh_series`（全日 120s/60s 值）寫入每日 result，之後唔使重跑先拎到峰值時刻
- **雙頻帶對照（L2）**：新增 0.1–0.5 Hz band（火山研究實錘 lock 帶），含第二頻帶 coherence + Hilbert 相位鎖定 `c=cos(Δφ)` kappa
- **判據並行**：`alert`（v1 分佈抬升型）與 `spike_alert`（v2 瞬時+事件對齊型）獨立輸出，歷史 schema 向後相容

驗證（本機重算）：8/25 M5.5 → `spike_alert=True`（z=89.7，峰值 06:59 UTC = 地震前 1 分鐘）；8/15 M4.9 → `spike_alert=True`（峰值 11:29 UTC = 地震前 1 分鐘）；8/01 無地震 → `spike_alert=False`（無誤報）。兩次事件峰值均出現於地震前 1 分鐘，pattern 一致。

## v3 升級（2026-11-13）

背景：2022-09-14 單控制日 baseline 已 stale。用 **10 個無震日**（TW 2026-10-02/08/14/19/25/31、11-05/07/10/12；CWA 有感目錄零事件；避開 10/22–23 颱風核心日）重校，pooled **23,510 窗 × 3 指標**，FP 驗證 0/10：

| 常數 | v2（2022-09-14） | v3（2026-10/11 重校） | 說明 |
|------|------|------|------|
| CONTROL_P99 | 0.205 | **0.145** | pooled p99；舊值偏高損靈敏度 |
| CONTROL_MEAN | 0.068 | **0.066** | 實質唔變 |
| ABS_SPIKE | 0.55 | **0.31** | pooled p99.9；舊值已跌入無震日背景（0.30–0.58）；0.31 保留 M3.x 同震偵測 |
| Z_SPIKE | 6.0 | 6.0（**裝飾性**） | 無震日 in-day z 39–83 恒過線（重尾），真 gate 係事件對齊 |

- **rolling 30-day 自校準**：每日由 `seismic_history.json` 用 median-of-daily-p99 自動追蹤背景漂移
  （10 無震日驗證 vs pooled 真值差 0.7%）；污染防護 = 排除事件日 + clamp [0.5×,2.5×]；
  乾淨日 <7 → fallback 固定常數。`baseline_source` 記錄於每日 result
- **v1 run 門檻 8→10**（`--run-threshold`）：無震日最長 run=5，margin 1.6×→2×
- **schema 延伸**：history record 新增 `band01_05_p99` / `band01_05_kappa_p99`
  （band B / kappa 嘅 rolling 等 30 日新數據累積後啟用；暫用固定參考線 0.133 / 0.183）

## ⚠️ 局限性（依然存在 — 2026-11-13 更新）

1. **M3.x 前兆測唔到**：2026-11-13 花蓮 M3.1/M3.8 backtest（pre-EQ 24h，三指標 × 兩頻帶）
   震前全 null，只有同震。能量差使然（M3.8 vs M7.4 ≈ 10^5×）—— **唔係框架失效，係量級
   未到門檻**。預期可偵測門檻 ~M4.5+（spike 驗證）/ M5.5+（前兆結構），待事件驗證
2. **觀測域限制**：本 pipeline 係**遠場 coherence**；近場鑽井孔隙水壓類前兆
   （如 T-18.1 分鐘 phase transition）屬另一觀測域，本 pipeline 未測試、null 唔構成反證
3. **單尖峰唔係證據**：背景重尾 —— 無震日日日都有 ≥0.30 尖峰（max 0.30–0.58）；
   **任何異常掃描必須配無震日對照**，否則會將背景波動誤讀成前兆（方法論鐵律）
4. **樣本量**：無震日 FP 驗證 n=10（颱風日排除）；新 baseline 未實測任何 M4.5+ 事件；
   M3.x 前兆 null 只有 n=2
5. **rolling 自校準係新嘢**：跨季節漂移追蹤行為未驗證；self-poisoning 防護（事件日
   排除 + clamp）已實作但長期表現待觀察
6. **Z_SPIKE 裝飾性**：v2 spike 偵測嘅實際 gate 一直係「事件對齊」而非 z 門檻（v4 可移除 z 路徑）
7. **數據收集進行中** —— 目前狀態：等 M4.5+ 事件做前瞻驗證（spike_alert + 新 baseline），
   同埋累積 band B / kappa 新 schema 數據。**呢條線而家係「watch and wait」模式**

## 使用

```bash
# 指定日期（前一日 24h）
python seismic_monitor.py --date 20260801

# 最近 24h（預設）
python seismic_monitor.py

# 自訂窗口 / 站
python seismic_monitor.py --hours 6 --stations NACB,YHNB,SSLB

# 自訂 fallback baseline / run 門檻
python seismic_monitor.py --baseline-p99 0.145 --run-threshold 10
```

## 輸出

- `seismic_result_YYYYMMDD.json` — 當日測量
- `seismic_history.json` — 累積時間序列

## 警報邏輯（v2 雙警報）

**`alert=True`（v1 分佈抬升型，對應火山 tremor 式持續抬升）**：
>baseline p99 比例 ≥5% **且** 最長連續 ≥10 窗口（v2 原值 8，v3 加固至 10；=10+ 分鐘持續抬升）。
baseline p99 v3 起用 rolling 30-day 自校準（fallback 0.145）。有 M≥5 地震時高 coherence
係預期行為；無地震時高 coherence 先係異常（待查）。

**`spike_alert=True`（v2 瞬時+對齊型，對應地震 spike）**：
全日 max_coh 為 spike（robust z ≥ 6〔裝飾性〕或 max_coh ≥ 0.31〔v3；v2 原值 0.55〕）**且** 與 USGS M≥4.5 事件 ±30 分鐘對齊。
v1 對瞬時 spike 免疫（8/25 M5.5 個案），v2 補返呢個盲區。

## 作者

tygtDc, Deep Research（D₁ 框架跨域研究）
