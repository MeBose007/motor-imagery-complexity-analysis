"""Statistics: Tables II-III (main text) and Tables S1-S2 (supplement)."""
import itertools
import numpy as np
import pandas as pd
from scipy import stats

from .aggregate import select

# ------------------------------------------------------------------ cohorts on PE_sub
def cohort_values(sub, group, condition):
    """PE_sub values (indexed by subject) for group in {minor, severe, LHA, RHA, healthy}."""
    if group in ("minor", "severe"):
        s = select(sub, condition, severity=group)
    elif group in ("LHA", "RHA"):
        s = select(sub, condition, lesion=group)
    else:
        s = select(sub, condition)
    return s.set_index("subject")["pe"]


def table2(sub):
    """Median and mean +/- SD of PE_sub (SD with ddof=0, as numpy's default)."""
    rows = []
    for grp, cond in [("minor", "L-MIT"), ("minor", "R-MIT"), ("severe", "L-MIT"),
                      ("severe", "R-MIT"), ("healthy", "EO-MI"), ("healthy", "EC-RS")]:
        v = cohort_values(sub, grp, cond).to_numpy()
        rows.append({"Group": grp.capitalize(), "Category": cond, "n": len(v),
                     "Median": np.median(v), "Mean": v.mean(), "SD": v.std()})
    return pd.DataFrame(rows)


TABLE3_PAIRS = [
    # (group1, condition1, group2, condition2, test)
    ("minor", "L-MIT", "severe", "L-MIT", "mannwhitney"),
    ("minor", "L-MIT", "severe", "R-MIT", "mannwhitney"),
    ("minor", "R-MIT", "severe", "L-MIT", "mannwhitney"),
    ("minor", "R-MIT", "severe", "R-MIT", "mannwhitney"),
    ("minor", "L-MIT", "healthy", "EO-MI", "mannwhitney"),
    ("minor", "L-MIT", "healthy", "EC-RS", "mannwhitney"),
    ("minor", "R-MIT", "healthy", "EO-MI", "mannwhitney"),
    ("minor", "R-MIT", "healthy", "EC-RS", "mannwhitney"),
    ("severe", "L-MIT", "healthy", "EO-MI", "mannwhitney"),
    ("severe", "L-MIT", "healthy", "EC-RS", "mannwhitney"),
    ("severe", "R-MIT", "healthy", "EO-MI", "mannwhitney"),
    ("severe", "R-MIT", "healthy", "EC-RS", "mannwhitney"),
    ("healthy", "EO-MI", "healthy", "EC-RS", "mannwhitney"),
    ("minor", "L-MIT", "minor", "R-MIT", "wilcoxon"),
    ("severe", "L-MIT", "severe", "R-MIT", "wilcoxon"),
]


def compare(sub, g1, c1, g2, c2, test):
    a, b = cohort_values(sub, g1, c1), cohort_values(sub, g2, c2)
    if test == "wilcoxon":               # paired by subject
        a, b = a.align(b, join="inner")
        res = stats.wilcoxon(a.to_numpy(), b.to_numpy(), alternative="two-sided")
    else:
        res = stats.mannwhitneyu(a.to_numpy(), b.to_numpy(), alternative="two-sided")
    return res.statistic, res.pvalue


def table3(sub, pairs=TABLE3_PAIRS):
    rows = []
    for g1, c1, g2, c2, test in pairs:
        stat, p = compare(sub, g1, c1, g2, c2, test)
        rows.append({"Group 1": g1.capitalize(), "Category 1": c1, "Group 2": g2.capitalize(),
                     "Category 2": c2, "Test": "Mann-Whitney U" if test == "mannwhitney" else "Wilcoxon signed-rank",
                     "Statistic": stat, "p-value": p})
    return pd.DataFrame(rows)


LESION_PAIRS = [
    ("LHA", "L-MIT", "LHA", "R-MIT", "wilcoxon"),
    ("RHA", "L-MIT", "RHA", "R-MIT", "wilcoxon"),
    ("LHA", "L-MIT", "RHA", "L-MIT", "mannwhitney"),
    ("LHA", "R-MIT", "RHA", "R-MIT", "mannwhitney"),
]


def format_p(p):
    return "<1e-6" if p < 1e-6 else f"{p:.2g}" if p < 0.01 else f"{p:.2f}"


# ------------------------------------------------------------------ Tables S1-S2
def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values (identical to statsmodels 'fdr_bh').
    Non-finite entries are left untouched; nothing is adjusted if <2 finite values."""
    p = np.asarray(pvals, dtype=float)
    out = p.copy()
    ok = np.isfinite(p)
    if ok.sum() < 2:
        return out
    q = p[ok]
    n = q.size
    order = np.argsort(q)
    ranked = q[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(n)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


COHORTS4 = ["minor", "severe", "EO", "EC"]
PAIR_LABELS = {("minor", "severe"): "M-S", ("minor", "EO"): "M-EO", ("minor", "EC"): "M-EC",
               ("severe", "EO"): "S-EO", ("severe", "EC"): "S-EC", ("EO", "EC"): "EO-EC"}


def trial_level_frame(df_published):
    """Per-trial PE in long format with a 4-cohort label. Stroke contributes all 40
    trials (left and right as recorded), healthy subjects their 20 blocks. Values are
    rounded to 6 decimals as in the authors' intermediate CSV."""
    d = df_published.copy()
    d["cohort"] = np.where(d["dataset"] == "stroke", d["severity"],
                           np.where(d["dataset"] == "eo_mi", "EO", "EC"))
    d["pe"] = d["pe"].round(6)
    return d[["cohort", "subject", "trial", "channel", "pe"]]


def _vectors(frame, channels, level):
    """One vector per cohort: all (subject, trial[, channel]) values when level='trial';
    one value per subject (mean over trials and channels) when level='subject'."""
    sel = frame[frame["channel"].isin(channels)]
    out = {}
    for c in COHORTS4:
        s = sel[sel["cohort"] == c]
        if level == "trial":
            out[c] = s["pe"].to_numpy()
        else:
            per_ch = s.groupby(["subject", "channel"])["pe"].mean()
            out[c] = per_ch.groupby("subject").mean().to_numpy()
    return out


def anova_posthoc(vectors):
    f, p = stats.f_oneway(*[vectors[c] for c in COHORTS4])
    res = {"F": f, "ANOVA_p": p}
    for (a, b), lab in PAIR_LABELS.items():
        res[f"{lab} Student"] = stats.ttest_ind(vectors[a], vectors[b], equal_var=True).pvalue
        res[f"{lab} Welch"] = stats.ttest_ind(vectors[a], vectors[b], equal_var=False).pvalue
    return res


def anova_tables(frame, canonical, regions, level="trial"):
    """Sensor-wise (Table S1) and region-wise (Table S2) one-way ANOVA with Student's
    and Welch's post-hoc t-tests; BH-FDR applied across sensors (or regions) separately
    for every column."""
    sensor = pd.DataFrame([{"Sensor": ch, **anova_posthoc(_vectors(frame, [ch], level))}
                           for ch in canonical])
    region = pd.DataFrame([{"Region": r, "Sensors": ", ".join(chs),
                            **anova_posthoc(_vectors(frame, chs, level))}
                           for r, chs in regions.items()])
    for tab in (sensor, region):
        for col in [c for c in tab.columns if c == "ANOVA_p" or c.endswith(("Student", "Welch"))]:
            tab[col + " (FDR)"] = bh_fdr(tab[col])
    return sensor, region
