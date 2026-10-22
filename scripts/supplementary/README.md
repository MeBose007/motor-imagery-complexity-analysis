# Supplementary material

| Output | Script | Notes |
|---|---|---|
| Fig. S1 | - | electrode montage diagram, no code |
| Fig. S2 | `figS2_tac.py` | `--include-white-noise` adds a white-noise reference (not in the published figure) |
| Fig. S3 | `figS3_alpha_band.py` | runs steps 1-3 with `--band alpha`; accepts `--mode` |
| Tables S1-S2 | `tables_S1_S2_anova.py` | needs `01_compute_pe.py --mode published` first |
| Table S3 | `table_S3_cbpt.py` | needs the same, plus `electrodes_tsv` in `config/paths.yaml`; the full sweep (N up to 200,000, k = 2-4) takes a while |
| Figs. S4-S6 | not yet included | resting-state PE of the Dortmund cohort |
| Figs. S7-S8 | `figS7_S8_psd.py` | Dortmund cohort, first 180 s |
| Table S4 | - | clinical metadata of the stroke dataset |

Tables S1-S3 follow Supplementary Sec. I D: PE values of all 4-s segments are pooled
across subjects within each cohort, so segments are treated as independent observations.
`tables_S1_S2_anova.py --level subject` gives the same analysis with one value per subject.
