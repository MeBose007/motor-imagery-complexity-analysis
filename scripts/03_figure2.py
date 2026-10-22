"""Step 3 - Fig. 2: PE_sub for (a) minor and (b) severe stroke against the two healthy
cohorts, and (c) L-MIT vs R-MIT in minor and severe stroke.
With --band alpha this produces Fig. S3."""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
import numpy as np
import matplotlib.pyplot as plt
from strokepe.cli import parser, setup, load_pe_table, save_figure
from strokepe.aggregate import pe_sub
from strokepe import stats, plotting

YLIM = {"broadband": (0.2, 0.8), "alpha": (0.2, 0.5)}
CONDS = ["L-MIT", "R-MIT", "EO-MI", "EC-RS"]


def panel_ab(ax, sub, group, ylim):
    groups = {c: stats.cohort_values(sub, group if c.endswith("MIT") else "healthy", c).to_numpy()
              for c in CONDS}
    plotting.violin_box(ax, groups, plotting.COLORS, ylim=ylim)
    pairs = []
    for i, j in [(0, 2), (1, 2), (2, 3), (0, 3), (1, 3)]:
        g1 = group if CONDS[i].endswith("MIT") else "healthy"
        _, p = stats.compare(sub, g1, CONDS[i], "healthy", CONDS[j], "mannwhitney")
        pairs.append((i, j, p))
    plotting.significance_bars(ax, pairs, groups)
    return groups


def panel_c(ax, sub):
    keys = [("minor", "L-MIT"), ("minor", "R-MIT"), ("severe", "L-MIT"), ("severe", "R-MIT")]
    groups = {f"{g}|{c}": stats.cohort_values(sub, g, c).to_numpy() for g, c in keys}
    colors = {k: plotting.COLORS[k.split("|")[1]] for k in groups}
    allv = np.concatenate(list(groups.values()))
    pad = (allv.max() - allv.min()) * 0.25
    plotting.violin_box(ax, groups, colors, ylim=(allv.min() - pad, allv.max() + pad * 2))
    pairs = [(0, 1, stats.compare(sub, "minor", "L-MIT", "minor", "R-MIT", "wilcoxon")[1]),
             (2, 3, stats.compare(sub, "severe", "L-MIT", "severe", "R-MIT", "wilcoxon")[1])]
    plotting.significance_bars(ax, pairs, groups, step=0.0)
    plotting.group_brackets(ax, [(0, 1, "minor"), (2, 3, "severe")])
    ax.set_ylabel("")


def main(argv=None):
    args = parser(__doc__).parse_args(argv)
    cfg, out = setup(args)
    sub = pe_sub(load_pe_table(cfg, args.mode, args.band))
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.5), gridspec_kw={"width_ratios": [1, 1, 0.9]})
    panel_ab(axes[0], sub, "minor", YLIM[args.band])
    panel_ab(axes[1], sub, "severe", YLIM[args.band])
    axes[1].set_ylabel("")
    panel_c(axes[2], sub)
    for ax, lab in zip(axes, "abc"):
        ax.set_title(f"({lab})", loc="left", fontsize=20)
    fig.tight_layout()
    name = "figure2" if args.band == "broadband" else "figureS3"
    save_figure(fig, out, name)


if __name__ == "__main__":
    main()
