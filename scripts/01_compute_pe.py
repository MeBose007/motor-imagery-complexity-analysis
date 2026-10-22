"""Step 1 - permutation entropy for every subject, condition and channel.

Writes results/<band>/<mode>/pe_table.csv (one row per subject x condition x channel
x segment), plus pe_sub.csv (PE_sub) and pe_sensor.csv (PE_sensor per group).
"""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
from strokepe.cli import parser, setup
from strokepe.compute import compute_pe_table
from strokepe.aggregate import pe_sub, pe_sensor
import pandas as pd


def sensor_table(df, canonical):
    cols = {}
    for name, cond, flt in [("minor L-MIT", "L-MIT", {"severity": "minor"}),
                            ("minor R-MIT", "R-MIT", {"severity": "minor"}),
                            ("severe L-MIT", "L-MIT", {"severity": "severe"}),
                            ("severe R-MIT", "R-MIT", {"severity": "severe"}),
                            ("LHA L-MIT", "L-MIT", {"lesion": "LHA"}),
                            ("LHA R-MIT", "R-MIT", {"lesion": "LHA"}),
                            ("RHA L-MIT", "L-MIT", {"lesion": "RHA"}),
                            ("RHA R-MIT", "R-MIT", {"lesion": "RHA"}),
                            ("EO-MI", "EO-MI", {}), ("EC-RS", "EC-RS", {})]:
        cols[name] = pe_sensor(df, cond, canonical, **flt)
    return pd.DataFrame(cols).rename_axis("channel")


def main(argv=None):
    p = parser(__doc__)
    p.add_argument("--jobs", type=int, default=1, help="parallel worker processes")
    args = p.parse_args(argv)
    cfg, out = setup(args)
    print(f"Computing PE  mode={args.mode}  band={args.band}")
    df = compute_pe_table(cfg, args.mode, args.band, jobs=args.jobs)
    df.to_csv(out / "pe_table.csv", index=False, float_format="%.10g")
    pe_sub(df).to_csv(out / "pe_sub.csv", index=False, float_format="%.10g")
    sensor_table(df, cfg["channels"]["canonical"]).to_csv(out / "pe_sensor.csv", float_format="%.10g")
    print(f"  wrote {out / 'pe_table.csv'} ({len(df)} rows), pe_sub.csv, pe_sensor.csv")


if __name__ == "__main__":
    main()
