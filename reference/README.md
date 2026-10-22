# Reference values

Numbers from the published paper and from the authors' analysis scripts, used by
`scripts/06_check_published.py` and the tests.

| File | Content | Source |
|---|---|---|
| `published_tables.yaml` | Table II (median, mean, SD) and Table III (p-values) | the paper |
| `published_subject_pe.csv` | PE_sub of the 18 subjects in each stroke group/condition, and the healthy values the authors used for the Table III tests | authors' statistics script |
| `published_sensor_pe.csv` | PE_sensor (25 channels) behind the topographic maps of Figs. 3(c-h) and 4(b-e) | authors' topoplot script |

Notes:

* In `published_subject_pe.csv` the stroke values reproduce Table II exactly. The
  healthy EO-MI/EC-RS values there come from a different run than Table II (their mean is
  0.6374 / 0.2965 instead of 0.6368 / 0.3035); the healthy means in Table II equal the
  means of the channel values in `published_sensor_pe.csv`.
* The values are listed in group order, not keyed by subject. For severe R-MIT, positions
  16-18 appear rotated relative to severe L-MIT (re-pairing them makes every L/R
  difference smaller than 0.005, like all other subjects). With that ordering the paired
  Wilcoxon test gives p = 0.18 (Table III); with the re-paired values p = 0.25. The
  pipeline pairs by subject ID, so expect 0.25 in `--mode published`.
