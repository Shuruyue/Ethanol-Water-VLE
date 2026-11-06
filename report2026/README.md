# report2026 — revised final report (reproducible)

Consulting report (XeLaTeX) on the non-ideal VLE of ethanol–water at 1 atm, its excess properties
and the energy demand of distillation. The report itself is written in Traditional Chinese, as
required by the course.
Its structure follows the course's Part I–V exactly (with a mapping to the grading weights); the
main text is 8 pages (excluding cover, table of contents and appendices A–G), within the required
6–10 pages excluding appendices.
**Every number** in the report is generated from `results/results.json` (via `tex/numbers.tex`), so
the text cannot drift from the code.

| Path | Contents |
|---|---|
| `calc/vle_core.py` | Vapour pressure (NIST Ambrose & Sprake Antoine for ethanol, IAPWS-95 for water), NRTL / Van Laar, Tsonopoulos virial correction, bubble / dew / azeotrope points, McCabe–Thiele and rating |
| `calc/run_analysis.py` | Parameter regression, experimental validation, Herington consistency test, excess properties, azeotrope vs pressure, pressure-swing distillation, column design, ethylene glycol ternary data → `results/` |
| `calc/make_figures.py` / `make_numbers.py` / `make_tables.py` | Vector PDF figures, LaTeX number macros, tables |
| `calc/test_core.py` | Consistency tests: Gibbs–Helmholtz analytic derivative, pure-component limits, azeotrope vs experiment |
| `data/` | Binary and ternary data of Kamihama et al., *J. Chem. Eng. Data* 57 (2012) 339 (NIST TRC ThermoML) |
| `VLE_Consulting_Report_2026.pdf` | Compiled report |

## Reproduce

```bash
pip install numpy scipy pandas matplotlib chemicals CoolProp pytest
cd calc
python run_analysis.py && python make_figures.py && python make_numbers.py && python make_tables.py
pytest -q test_core.py
cd ../tex && xelatex report.tex && xelatex report.tex   # needs ctex (TeX Live: texlive-lang-chinese)
```

The group members on the cover are set by the `\MemberList` macro near the top of `tex/report.tex`.
