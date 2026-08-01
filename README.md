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

## 警報邏輯

`alert=True` 當：>baseline p99 比例 ≥5% **且** 最長連續 ≥8 窗口（=8+ 分鐘持續抬升）。
有 M≥5 地震時高 coherence 係預期行為；無地震時高 coherence 先係異常（待查）。

## 作者

tygtDc, Deep Research（D₁ 框架跨域研究）
