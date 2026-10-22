# StrokePE

[![DOI](https://zenodo.org/badge/1429544572.svg)](https://doi.org/10.5281/zenodo.23559636)

Code for

> A. Bose, P. Gupta, K. Vemuri, D. Das, D. Ghosh, S. K. Dana, S. Saha, and C. Hens,
> "Ordinal pattern of brain electrical activity as a marker of stroke-induced alterations
> in motor imagery task," *Chaos* **36**, 063101 (2026). https://doi.org/10.1063/5.0314904

Permutation entropy (PE) of EEG from 36 acute stroke patients performing left- and
right-hand motor imagery (L-MIT, R-MIT) is compared with healthy eye-open motor imagery
(EO-MI) and healthy eye-closed rest (EC-RS).

## Two modes

The paper's Methods text and the computation behind the published figures differ in one
step, so the pipeline offers both:

| | `--mode paper` (default) | `--mode published` |
|---|---|---|
| Follows | Sec. II of the paper | the analysis that produced the published figures and tables |
| Stroke | the twenty 4-s motor-imagery segments per hand are concatenated into one 80-s series; PE once per channel | PE of each 4-s segment, averaged over the 20 segments |
| Healthy | PE of the first 80 s of each recording | first 80 s split into 20 blocks of 4 s; PE per block, averaged |
| Fig. 3(a, b) | one value per subject in each box | every subject's 4-s segments pooled in each box |

Everything else is shared: subjects and groups, the 25 channels, mu = 5, tau = 1,
normalisation by log(mu!), PE_sub (mean over channels) and PE_sensor (mean over subjects),
the statistical tests and the figures. Use `--mode published` to regenerate the published
numbers; `scripts/06_check_published.py` compares them with the paper.

## Quick start

```bash
pip install -r requirements.txt
cp config/paths.example.yaml config/paths.yaml     # then edit the data folders
python scripts/run_all.py                          # paper method
python scripts/run_all.py --mode published         # published numbers + check report
```

Outputs go to `results/<band>/<mode>/`. `--jobs N` runs subjects in parallel.

## What reproduces what

| Output | Command | Files |
|---|---|---|
| PE per subject, condition, channel | `scripts/01_compute_pe.py` | `pe_table.csv`, `pe_sub.csv`, `pe_sensor.csv` |
| Tables II and III | `scripts/02_tables.py` | `table2.csv`, `table3.csv`, `tables.md` |
| Fig. 2 | `scripts/03_figure2.py` | `figure2.png/.pdf` |
| Fig. 3 | `scripts/04_figure3.py` | `figure3ab`, `figure3c-h` (+ value ranges) |
| Fig. 4 | `scripts/05_figure4.py` | `figure4a`, `figure4b-e`, `fig4_lesion_tests.csv` |
| Check against the paper | `scripts/06_check_published.py --mode published` | `check_report.md` |
| Fig. S2 | `scripts/supplementary/figS2_tac.py` | `results/supplementary/figureS2*` |
| Fig. S3 | `scripts/supplementary/figS3_alpha_band.py` | `results/alpha/<mode>/figureS3*` |
| Tables S1-S2 | `scripts/supplementary/tables_S1_S2_anova.py` | `results/broadband/supplementary/` |
| Table S3 | `scripts/supplementary/table_S3_cbpt.py` | `results/broadband/supplementary/` |
| Figs. S7-S8 | `scripts/supplementary/figS7_S8_psd.py` | `results/supplementary/` |

See `scripts/supplementary/README.md` for details. Figs. S4-S6 (resting-state PE of the
Dortmund cohort) are not included yet.

## Data

All datasets are public; see `data/README.md` for sources and expected file names, and
`preprocessing/README.md` for the EEGLAB preprocessing.

## Method as implemented

* **Stroke** (500 Hz): each recording holds 40 trials of 8 s (2 s cue, 4 s imagery, 2 s
  rest) alternating left/right, starting with left. Samples 1000-2999 of each trial are the
  motor-imagery segment. 25 channels common to all datasets (Fig. S1).
* **EO-MI** (PhysioNet run R08, 160 Hz) and **EC-RS** (OpenNeuro ds003775, 1024 Hz, cleaned
  epochs joined end to end): first 80 s. No resampling.
* **PE**: ordinal patterns of length mu = 5, delay tau = 1; ties broken by temporal order;
  Shannon entropy of the pattern distribution divided by log(5!).
* **PE_sub** = mean over the 25 channels; **PE_sensor** = mean over subjects per channel;
  topographic maps are min-max normalised per condition.
* **Statistics** (two-sided, alpha = 0.05): Mann-Whitney U between groups, Wilcoxon
  signed-rank between L-MIT and R-MIT of the same subjects (paired by subject).
  Table II reports mean and SD with ddof = 0.
* **Fig. S3**: 8-13 Hz, 4th-order Butterworth applied forwards and backwards; stroke and
  EO-MI filtered as continuous recordings, EC-RS epoch by epoch.

## Notes on the published version

*Draft for the authors to review before release.*

1. **Methods text vs computation.** Sec. II describes PE on 80-s series; the published
   figures and tables used PE per 4-s segment, averaged (see "Two modes").
2. **Table III, severe L-MIT vs R-MIT.** The published p = 0.18 was computed from
   per-subject values in which the last three severe R-MIT entries were out of order
   relative to L-MIT. Pairing by subject, as the pipeline does, gives p of about 0.25; the
   conclusion (not significant) is unchanged. See `reference/README.md`.
3. **Severity groups.** The minor/severe lists in `config/subjects.yaml` reproduce Table I.
   Six patients in the severe group have NIHSS <= 4 (Table S4 IDs 12, 18, 20, 22, 23, 29),
   so the grouping did not follow NIHSS alone.
4. **Fig. 4.** The caption ((b, c) = left hemisphere affected) matches the computed maps;
   the text of Sec. III B refers to (b, c) as right-lesioned.
5. **Fig. 3 caption ranges.** The computed ranges are EC-RS [0.28, 0.32] and minor L-MIT
   [0.35, 0.37] (printed as [0.23, 0.32] and [0.35, 0.32]).
6. **Tables S1-S3** pool the PE of all 4-s segments of all subjects in a cohort and treat
   them as independent observations (stroke: 40 trials per subject with left and right
   together; healthy: 20 blocks). `tables_S1_S2_anova.py --level subject` gives the same
   analysis with one value per subject.
7. **Table S3.** Cluster-forming threshold alpha = 0.05 two-sided. Permutation p-values
   cannot be smaller than 1/(N+1); entries printed as 0 correspond to p < 5e-6
   (N = 200,000).

## Tests

```bash
python -m pytest tests
```

The tests check the PE and a_c implementations against the authors' original functions,
recompute Table II (stroke rows) and Table III from the published per-subject values, and
run every script on a small synthetic dataset in both modes.

## License

MIT (see `LICENSE`). Please cite the article above (`CITATION.cff`).
