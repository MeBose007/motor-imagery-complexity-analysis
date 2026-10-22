"""Small helpers shared by the scripts in scripts/."""
import argparse
from pathlib import Path
import pandas as pd

from .config import load_config, results_dir

MODE_HELP = ("'paper' (default): PE as described in the paper (80-s series). "
             "'published': PE per 4-s segment, averaged - the computation behind the "
             "published figures and tables.")


def parser(description, mode=True, band=True):
    p = argparse.ArgumentParser(description=description)
    if mode:
        p.add_argument("--mode", choices=["paper", "published"], default="paper", help=MODE_HELP)
    if band:
        p.add_argument("--band", choices=["broadband", "alpha"], default="broadband",
                       help="broadband (main text) or alpha = 8-13 Hz (Fig. S3)")
    p.add_argument("--paths", default=None, help="paths YAML (default: config/paths.yaml)")
    return p


def setup(args):
    cfg = load_config(args.paths)
    out = results_dir(cfg, getattr(args, "band", "broadband"), getattr(args, "mode", None))
    return cfg, out


def pe_table_path(cfg, mode, band):
    return results_dir(cfg, band, mode) / "pe_table.csv"


def load_pe_table(cfg, mode, band):
    path = pe_table_path(cfg, mode, band)
    if not path.exists():
        raise SystemExit(f"{path} not found - run: python scripts/01_compute_pe.py "
                         f"--mode {mode} --band {band}")
    return pd.read_csv(path, dtype={"subject": str})


def markdown(df, floatfmt="{:.4f}"):
    """Minimal DataFrame -> Markdown table (no extra dependency)."""
    def fmt(v):
        if isinstance(v, float):
            return floatfmt.format(v)
        return str(v)
    head = "| " + " | ".join(map(str, df.columns)) + " |"
    sep = "|" + "|".join("---" for _ in df.columns) + "|"
    body = ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([head, sep, *body])


def save_figure(fig, out_dir, name):
    out_dir = Path(out_dir)
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"{name}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    print(f"  saved {out_dir / name}.png/.pdf")
