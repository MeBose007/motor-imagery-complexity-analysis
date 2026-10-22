"""Step 2 - Table II (median, mean +/- SD of PE_sub) and Table III (Mann-Whitney U and
Wilcoxon signed-rank tests), plus the lesion-side comparisons behind Fig. 4(a)."""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
from strokepe.cli import parser, setup, load_pe_table, markdown
from strokepe.aggregate import pe_sub
from strokepe import stats


def main(argv=None):
    args = parser(__doc__).parse_args(argv)
    cfg, out = setup(args)
    sub = pe_sub(load_pe_table(cfg, args.mode, args.band))
    t2 = stats.table2(sub)
    t3 = stats.table3(sub)
    t4 = stats.table3(sub, stats.LESION_PAIRS)
    t2.to_csv(out / "table2.csv", index=False)
    t3.to_csv(out / "table3.csv", index=False)
    t4.to_csv(out / "fig4_lesion_tests.csv", index=False)
    show3 = t3.assign(**{"p-value": t3["p-value"].map(stats.format_p)}).drop(columns="Statistic")
    show4 = t4.assign(**{"p-value": t4["p-value"].map(stats.format_p)}).drop(columns="Statistic")
    text = (f"# Tables (mode={args.mode}, band={args.band})\n\n## Table II\n\n{markdown(t2)}\n\n"
            f"## Table III\n\n{markdown(show3)}\n\n## Lesion side (Fig. 4a)\n\n{markdown(show4)}\n")
    (out / "tables.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
