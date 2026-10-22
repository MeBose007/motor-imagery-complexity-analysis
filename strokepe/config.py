"""Load the YAML configuration in config/."""
import os
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
# STROKEPE_CONFIG_DIR lets the tests (or a user) point to another set of YAML files.
CONFIG_DIR = Path(os.environ.get("STROKEPE_CONFIG_DIR", ROOT / "config"))


def _load(path):
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def load_config(paths_file=None):
    """Return a dict with keys: subjects, channels, params, paths."""
    config_dir = Path(os.environ.get("STROKEPE_CONFIG_DIR", CONFIG_DIR))
    cfg = {name: _load(config_dir / f"{name}.yaml") for name in ("subjects", "channels", "params")}
    pf = Path(paths_file) if paths_file else config_dir / "paths.yaml"
    if not pf.exists():
        raise FileNotFoundError(
            f"{pf} not found. Copy config/paths.example.yaml to config/paths.yaml "
            "and point it to your data folders.")
    raw = _load(pf) or {}
    base = pf.resolve().parent

    def resolve(value):
        p = Path(value).expanduser()
        return p if p.is_absolute() else (base / p).resolve()

    cfg["paths"] = {k: resolve(v) for k, v in raw.items() if v}
    return cfg


def results_dir(cfg, band="broadband", mode=None, create=True):
    """results/<band>/<mode>/ (or results/<band>/ when mode is None)."""
    d = Path(cfg["paths"].get("results_dir", ROOT / "results")) / band
    if mode:
        d = d / mode
    if create:
        d.mkdir(parents=True, exist_ok=True)
    return d


def subject_file(cfg, dataset, sid, **fmt):
    """Full path of a subject's input file."""
    sub = cfg["subjects"][dataset]
    folder = cfg["paths"][f"{dataset}_dir"]
    return Path(folder) / sub["file_pattern"].format(id=sid, **fmt)


def severity_of(cfg, sid):
    for grp, ids in cfg["subjects"]["stroke"]["severity"].items():
        if sid in ids:
            return grp
    return "-"


def lesion_of(cfg, sid):
    for grp, ids in cfg["subjects"]["stroke"]["lesion_side"].items():
        if sid in ids:
            return grp
    return "-"
