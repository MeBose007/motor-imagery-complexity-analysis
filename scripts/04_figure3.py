"""Step 4 - Fig. 3: (a, b) channel-wise PE for minor and severe stroke with the healthy
cohorts; (c-h) topographic maps of PE_sensor, each min-max normalised.
In 'paper' mode each box holds one value per subject; in 'published' mode it pools every
subject's 4-s segments, as in the published figure."""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
import matplotlib.pyplot as plt
import pandas as pd
from strokepe.cli import parser, setup, load_pe_table, save_figure
from strokepe.aggregate import subject_channel, pe_sensor
from strokepe import plotting

CONDS = ["L-MIT", "R-MIT", "EO-MI", "EC-RS"]
MAPS = [("(c) minor L-MIT", "L-MIT", {"severity": "minor"}),
        ("(d) minor R-MIT", "R-MIT", {"severity": "minor"}),
        ("(e) severe L-MIT", "L-MIT", {"severity": "severe"}),
        ("(f) severe R-MIT", "R-MIT", {"severity": "severe"}),
        ("(g) EO-MI", "EO-MI", {}), ("(h) EC-RS", "EC-RS", {})]


def box_data(df, severity, canonical):
    data = {}
    for cond in CONDS:
        d = df[df["condition"] == cond]
        if cond.endswith("MIT"):
            d = d[d["severity"] == severity]
        for ch in canonical:
            data[(cond, ch)] = d.loc[d["channel"] == ch, "pe"].to_numpy()
    return data


def draw_maps(df, canonical, maps, out, name, ncols=4):
    info = plotting.montage_info(canonical)
    nrows = -(-len(maps) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols + 1, 3.4 * nrows))
    ranges = []
    im = None
    for ax, (title, cond, flt) in zip(axes.flat, maps):
        v = pe_sensor(df, cond, canonical, **flt)
        im = plotting.topoplot(ax, v.to_numpy(), info)
        ax.set_title(title, fontsize=12)
        ranges.append({"map": title, "min": v.min(), "max": v.max()})
    for ax in list(axes.flat)[len(maps):]:
        ax.axis("off")
    cb = fig.colorbar(im, ax=list(axes.flat), shrink=0.6, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
    cb.set_label("normalised PE")
    save_figure(fig, out, name)
    pd.DataFrame(ranges).to_csv(out / f"{name}_ranges.csv", index=False, float_format="%.4f")
    return ranges


def main(argv=None):
    args = parser(__doc__).parse_args(argv)
    cfg, out = setup(args)
    canonical = cfg["channels"]["canonical"]
    df = load_pe_table(cfg, args.mode, args.band)
    if args.mode == "paper":
        df = subject_channel(df)          # one value per subject and channel
    fig, axes = plt.subplots(2, 1, figsize=(20, 10), sharex=True)
    for ax, sev, lab in zip(axes, ["minor", "severe"], "ab"):
        plotting.channel_boxplots(ax, box_data(df, sev, canonical), canonical, CONDS)
        ax.set_title(f"({lab}) {sev} stroke", loc="left", fontsize=16)
    fig.tight_layout()
    save_figure(fig, out, "figure3ab")
    ranges = draw_maps(df, canonical, MAPS, out, "figure3c-h", ncols=3)
    for r in ranges:
        print(f"  {r['map']:18s} range [{r['min']:.4f}, {r['max']:.4f}]")


if __name__ == "__main__":
    main()
