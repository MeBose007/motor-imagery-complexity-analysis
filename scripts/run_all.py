"""Run the main-text pipeline: PE, Tables II-III, Figs. 2-4 (and, in published mode,
the comparison with the published numbers).

    python scripts/run_all.py                    # paper method (default)
    python scripts/run_all.py --mode published   # regenerates the published numbers
    python scripts/run_all.py --band alpha       # 8-13 Hz analysis (Fig. S3)
"""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
import importlib.util
import pathlib
from strokepe.cli import parser

HERE = pathlib.Path(__file__).resolve().parent


def run(script, argv):
    spec = importlib.util.spec_from_file_location(script, HERE / script)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print(f"\n=== {script} {' '.join(argv)}")
    mod.main(argv)


def main():
    p = parser(__doc__)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--skip-compute", action="store_true", help="reuse an existing pe_table.csv")
    a = p.parse_args()
    common = ["--mode", a.mode, "--band", a.band] + (["--paths", a.paths] if a.paths else [])
    if not a.skip_compute:
        run("01_compute_pe.py", common + ["--jobs", str(a.jobs)])
    run("02_tables.py", common)
    run("03_figure2.py", common)
    if a.band == "broadband":
        run("04_figure3.py", common)
        run("05_figure4.py", common)
        if a.mode == "published":
            run("06_check_published.py", ["--mode", a.mode] + (["--paths", a.paths] if a.paths else []))


if __name__ == "__main__":
    main()
