"""Fig. S3 - the analysis of Fig. 2 repeated on 8-13 Hz band-passed signals
(4th-order Butterworth, zero-phase; stroke and EO-MI filtered as continuous recordings,
EC-RS epoch by epoch, as in the authors' filtering scripts)."""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
import importlib.util
import pathlib
from strokepe.cli import parser

SCRIPTS = pathlib.Path(__file__).resolve().parents[1]


def run(script, argv):
    spec = importlib.util.spec_from_file_location(script, SCRIPTS / script)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.main(argv)


def main():
    p = parser(__doc__, band=False)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    common = ["--mode", a.mode, "--band", "alpha"] + (["--paths", a.paths] if a.paths else [])
    run("01_compute_pe.py", common + ["--jobs", str(a.jobs)])
    run("02_tables.py", common)
    run("03_figure2.py", common)


if __name__ == "__main__":
    main()
