"""Step 6 - compare the pipeline output with the published numbers (use with
--mode published). Writes check_report.md next to the results."""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
import numpy as np
import pandas as pd
import yaml
from strokepe.cli import parser, setup, load_pe_table
from strokepe.config import ROOT
from strokepe.aggregate import pe_sub
from strokepe import stats

REF = ROOT / "reference"
# Published values that the pipeline is not expected to reproduce, with the reason.
KNOWN = {("table3", "severe L-MIT vs severe R-MIT"):
         "published p used values whose last three severe R-MIT entries were out of order "
         "(see reference/README.md); pairing by subject gives about 0.25"}


def check(rows, kind, name, ours, published, tol):
    ok = abs(ours - published) <= tol
    status = "PASS" if ok else ("EXPECTED DIFFERENCE" if (kind, name) in KNOWN else "FAIL")
    rows.append({"check": kind, "item": name, "pipeline": ours, "published": published, "status": status,
                 "note": KNOWN.get((kind, name), "")})


def main(argv=None):
    args = parser(__doc__, band=False).parse_args(argv)
    args.band = "broadband"
    cfg, out = setup(args)
    if args.mode != "published":
        print("Note: the published numbers come from --mode published; differences are expected in paper mode.")
    df = load_pe_table(cfg, args.mode, "broadband")
    sub = pe_sub(df)
    ref = yaml.safe_load(open(REF / "published_tables.yaml"))
    rows = []

    t2 = stats.table2(sub)
    for _, r in t2.iterrows():
        key = f"{r['Group'].lower()} {r['Category']}"
        med, mean, sd = ref["table2"][key]
        check(rows, "table2", f"{key} median", r["Median"], med, 5e-5)
        check(rows, "table2", f"{key} mean", r["Mean"], mean, 5e-5)
        check(rows, "table2", f"{key} SD", r["SD"], sd, 5e-5)

    for g1, c1, g2, c2, p_pub in ref["table3"]:
        test = "wilcoxon" if g1 == g2 and g1 != "healthy" else "mannwhitney"
        _, p = stats.compare(sub, g1, c1, g2, c2, test)
        name = f"{g1} {c1} vs {g2} {c2}"
        if p_pub <= 1e-6:
            check(rows, "table3", name, float(p < 1e-6), 1.0, 0)
            rows[-1].update(pipeline=p, published="<1e-6")
        else:
            check(rows, "table3", name, round(p, 2), p_pub, 0.005)

    subj = pd.read_csv(REF / "published_subject_pe.csv")
    for (g, c), ref_vals in subj[subj["group"] != "healthy"].groupby(["group", "condition"]):
        ours = np.sort(stats.cohort_values(sub, g, c).to_numpy())
        diff = np.max(np.abs(ours - np.sort(ref_vals["pe_sub"].to_numpy()))) if len(ours) == len(ref_vals) else np.inf
        check(rows, "PE_sub (18 subjects)", f"{g} {c} max |diff|", diff, 0.0, 1e-6)

    sensor_ref = pd.read_csv(REF / "published_sensor_pe.csv", index_col="channel")
    sensor = pd.read_csv(out / "pe_sensor.csv", index_col="channel")
    for col in sensor_ref.columns:
        diff = float(np.max(np.abs(sensor[col] - sensor_ref[col])))
        check(rows, "PE_sensor (25 channels)", f"{col} max |diff|", diff, 0.0, 1e-6)

    rep = pd.DataFrame(rows)
    lines = [f"# Check against published values (mode={args.mode})", "",
             "| check | item | pipeline | published | status |", "|---|---|---|---|---|"]
    for r in rep.itertuples():
        pv = f"{r.pipeline:.4g}" if isinstance(r.pipeline, float) else r.pipeline
        lines.append(f"| {r.check} | {r.item} | {pv} | {r.published} | {r.status} |")
    notes = sorted(set(n for n in rep["note"] if n))
    if notes:
        lines += ["", "Notes:", *[f"* {n}" for n in notes]]
    (out / "check_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print("\nSummary:", rep["status"].value_counts().to_dict())


if __name__ == "__main__":
    main()
