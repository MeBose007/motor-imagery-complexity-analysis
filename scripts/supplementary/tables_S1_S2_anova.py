"""Tables S1-S2 - sensor-wise and region-wise one-way ANOVA across the four cohorts,
with Student's and Welch's post-hoc t-tests and Benjamini-Hochberg FDR correction.

As described in Supplementary Sec. I D1, the default (--level trial) pools the PE of every
4-s segment of every subject in a cohort (stroke: all 40 trials; healthy: 20 blocks), so
segments are treated as independent observations. --level subject uses one value per
subject instead. Needs the per-segment PE table: run 01_compute_pe.py --mode published."""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
from strokepe.cli import parser, setup, load_pe_table, markdown
from strokepe import stats

LABELS = list(stats.PAIR_LABELS.values())


def short(tab, key, test):
    cols = [key, "ANOVA_p (FDR)"] + [f"{l} {test} (FDR)" for l in LABELS]
    t = tab[cols].copy()
    t.columns = [key, "ANOVA"] + LABELS
    for c in t.columns[1:]:
        t[c] = t[c].map(lambda p: f"{p:.2g}" + ("*" if p < 0.05 else ""))
    return t


def main(argv=None):
    p = parser(__doc__, mode=False)
    p.add_argument("--level", choices=["trial", "subject"], default="trial")
    args = p.parse_args(argv)
    cfg, out = setup(args)
    out = out / "supplementary"
    out.mkdir(exist_ok=True)
    frame = stats.trial_level_frame(load_pe_table(cfg, "published", args.band))
    sensor, region = stats.anova_tables(frame, cfg["channels"]["canonical"], cfg["channels"]["regions"], args.level)
    tag = f"{args.level}-level"
    sensor.to_csv(out / f"tableS1_sensor_anova_{tag}.csv", index=False)
    region.to_csv(out / f"tableS2_region_anova_{tag}.csv", index=False)
    md = [f"# Tables S1-S2 ({tag}, band={args.band}); FDR-corrected p, * = p_FDR < 0.05", ""]
    for test in ("Student", "Welch"):
        md += [f"## Table S1 - sensors, {test}'s t", "", markdown(short(sensor, "Sensor", test)), ""]
    md += ["## Table S2 - regions, Student's t", "", markdown(short(region, "Region", "Student")), ""]
    (out / f"tablesS1_S2_{tag}.md").write_text("\n".join(md), encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
