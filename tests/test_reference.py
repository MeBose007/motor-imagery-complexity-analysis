"""Table II (stroke rows) and Table III recomputed from the per-subject values in the
authors' statistics script (reference/published_subject_pe.csv)."""
import numpy as np
import pandas as pd
import yaml
from strokepe import stats
from strokepe.config import ROOT


def reference_sub():
    ref = pd.read_csv(ROOT / "reference/published_subject_pe.csv")
    rows = []
    for r in ref.itertuples():
        healthy = r.group == "healthy"
        rows.append({"dataset": "eo_mi" if r.condition == "EO-MI" else ("ec_rs" if healthy else "stroke"),
                     "subject": f"{r.condition if healthy else r.group}-{r.position:02d}",
                     "severity": "-" if healthy else r.group, "lesion": "-",
                     "condition": r.condition, "pe": r.pe_sub})
    return pd.DataFrame(rows)


def test_table2_stroke_rows():
    t2 = stats.table2(reference_sub())
    pub = yaml.safe_load(open(ROOT / "reference/published_tables.yaml"))["table2"]
    for _, r in t2.iterrows():
        if r["Group"] == "Healthy":
            continue
        med, mean, sd = pub[f"{r['Group'].lower()} {r['Category']}"]
        assert round(r["Median"], 4) == med and round(r["Mean"], 4) == mean and round(r["SD"], 4) == sd


def test_table3():
    t3 = stats.table3(reference_sub())
    pub = yaml.safe_load(open(ROOT / "reference/published_tables.yaml"))["table3"]
    for (_, r), (*_, p_pub) in zip(t3.iterrows(), pub):
        if p_pub <= 1e-6:
            assert r["p-value"] < 1e-6
        else:
            assert round(r["p-value"], 2) == p_pub


def test_bh_fdr():
    p = np.array([0.01, 0.04, 0.03, 0.005, np.nan])
    expected = np.array([0.02, 0.04, 0.04, 0.02, np.nan])
    np.testing.assert_allclose(stats.bh_fdr(p), expected, equal_nan=True)
