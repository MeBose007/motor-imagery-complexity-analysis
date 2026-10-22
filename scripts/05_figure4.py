"""Step 5 - Fig. 4: (a) PE_sub for L-MIT and R-MIT by lesion side (LHA / RHA = left /
right hemisphere affected); (b-e) topographic maps of PE_sensor."""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
import importlib.util
import pathlib
import matplotlib.pyplot as plt
from strokepe.cli import parser, setup, load_pe_table, save_figure
from strokepe.aggregate import pe_sub
from strokepe import stats, plotting

_spec = importlib.util.spec_from_file_location("fig3", pathlib.Path(__file__).with_name("04_figure3.py"))
fig3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fig3)

MAPS = [("(b) LHA L-MIT", "L-MIT", {"lesion": "LHA"}), ("(c) LHA R-MIT", "R-MIT", {"lesion": "LHA"}),
        ("(d) RHA L-MIT", "L-MIT", {"lesion": "RHA"}), ("(e) RHA R-MIT", "R-MIT", {"lesion": "RHA"})]


def main(argv=None):
    args = parser(__doc__).parse_args(argv)
    cfg, out = setup(args)
    df = load_pe_table(cfg, args.mode, args.band)
    sub = pe_sub(df)
    keys = [("LHA", "L-MIT"), ("LHA", "R-MIT"), ("RHA", "L-MIT"), ("RHA", "R-MIT")]
    groups = {f"{g}|{c}": stats.cohort_values(sub, g, c).to_numpy() for g, c in keys}
    colors = {k: plotting.COLORS[k.split("|")[1]] for k in groups}
    fig, ax = plt.subplots(figsize=(6, 6))
    allv = [v for vals in groups.values() for v in vals]
    pad = (max(allv) - min(allv)) * 0.25
    plotting.violin_box(ax, groups, colors, ylim=(min(allv) - pad, max(allv) + pad))
    plotting.group_brackets(ax, [(0, 1, "LHA"), (2, 3, "RHA")])
    ax.set_title("(a)", loc="left", fontsize=20)
    save_figure(fig, out, "figure4a")
    ranges = fig3.draw_maps(df, cfg["channels"]["canonical"], MAPS, out, "figure4b-e", ncols=2)
    for r in ranges:
        print(f"  {r['map']:16s} range [{r['min']:.4f}, {r['max']:.4f}]")


if __name__ == "__main__":
    main()
