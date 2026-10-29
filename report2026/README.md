# report2026 — 期末報告修正版（可重現）

乙醇–水 1 atm 非理想 VLE、過剩性質與蒸餾能耗的顧問報告（XeLaTeX）。
結構完全依課程規定的 Part I–V（含評分比重對照），正文 8 頁（不含封面、目錄、附錄 A–G），符合「6–10 頁（不含附錄）」。
報告中**每一個數字**都由 `results/results.json` 自動產生（`tex/numbers.tex`），避免文字與程式不一致。

| 路徑 | 內容 |
|---|---|
| `calc/vle_core.py` | 蒸氣壓（NIST Ambrose & Sprake 乙醇 Antoine、IAPWS-95 水）、NRTL／Van Laar、Tsonopoulos 維里修正、泡點／露點／共沸點、McCabe–Thiele 與 rating |
| `calc/run_analysis.py` | 參數迴歸、實驗驗證、Herington 一致性、過剩性質、共沸點 vs 壓力、PSD、塔設計、乙二醇三元數據分析 → `results/` |
| `calc/make_figures.py` / `make_numbers.py` / `make_tables.py` | 向量 PDF 圖、LaTeX 數字巨集、表格 |
| `calc/test_core.py` | Gibbs–Helmholtz 解析式、純物質極限、共沸點 vs 實驗等一致性測試 |
| `data/` | Kamihama et al., *J. Chem. Eng. Data* 57 (2012) 339 的二元與三元數據（NIST TRC ThermoML） |
| `VLE_Consulting_Report_2026.pdf` | 編譯好的報告 |

## 重現

```bash
pip install numpy scipy pandas matplotlib chemicals CoolProp pytest
cd calc
python run_analysis.py && python make_figures.py && python make_numbers.py && python make_tables.py
pytest -q test_core.py
cd ../tex && xelatex report.tex && xelatex report.tex   # 需 ctex（TeX Live: texlive-lang-chinese）
```

封面組員資訊在 `tex/report.tex` 開頭的 `\Members` 巨集。
