"""Permutation entropy for every subject, condition, channel (and segment).

Two modes:

* ``paper``      - the method as written in Sec. II: the twenty 4-s motor-imagery
                   segments per hand are concatenated into one 80-s series and PE is
                   computed once per channel; for healthy subjects PE is computed on the
                   first 80 s of each recording.
* ``published``  - the computation that produced the published figures and tables:
                   PE on each 4-s segment (20 per hand for stroke; 20 equal blocks of
                   the first 80 s for healthy subjects), averaged afterwards.
"""
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd

from . import io, segment, filtering
from .config import subject_file, severity_of, lesion_of
from .pe import permutation_entropy

COLUMNS = ["dataset", "subject", "severity", "lesion", "condition", "channel",
           "segment", "trial", "n_samples", "pe"]


def _bandpass(cfg, band, data, fs, axis=-1):
    if band == "broadband":
        return data
    lo, hi = cfg["params"]["bands"][band]
    return filtering.bandpass(data, lo, hi, fs, order=cfg["params"]["bands"]["filter_order"], axis=axis)


def _rows(base, condition, pieces, canonical, mu, tau):
    """pieces: list of (segment_index, trial_number, array channels x samples)."""
    rows = []
    for seg_idx, trial, arr in pieces:
        for ci, ch in enumerate(canonical):
            rows.append({**base, "condition": condition, "channel": ch, "segment": seg_idx,
                         "trial": trial, "n_samples": arr.shape[1],
                         "pe": permutation_entropy(arr[ci], mu, tau)})
    return rows


def _stroke(cfg, sid, mode, band):
    p, chans = cfg["params"], cfg["channels"]
    canonical = chans["canonical"]
    data = io.load_stroke_csv(subject_file(cfg, "stroke", sid), chans["stroke_drop_columns"], canonical)
    data = _bandpass(cfg, band, data, p["stroke"]["fs"])   # whole recording, as in the authors' filtering
    segs = segment.stroke_mi_segments(data, p["stroke"]["trial_samples"], p["stroke"]["mi_start"],
                                      p["stroke"]["mi_stop"], p["stroke"]["trials_expected"])
    base = {"dataset": "stroke", "subject": sid, "severity": severity_of(cfg, sid),
            "lesion": lesion_of(cfg, sid)}
    rows = []
    for cond in ("L-MIT", "R-MIT"):
        hand = [(t, s) for t, c, s in segs if c == cond]
        if mode == "paper":
            pieces = [(0, 0, np.concatenate([s for _, s in hand], axis=1))]
        else:
            pieces = [(k, t + 1, s) for k, (t, s) in enumerate(hand)]
        rows += _rows(base, cond, pieces, canonical, p["pe"]["mu"], p["pe"]["tau"])
    return rows


def _healthy_pieces(cfg, data, mode):
    if mode == "paper":
        return [(0, 0, data)]
    blocks = segment.equal_blocks(data, cfg["params"]["published_mode"]["n_blocks"])
    return [(k, k + 1, b) for k, b in enumerate(blocks)]


def _eo_mi(cfg, sid, mode, band):
    p, chans = cfg["params"], cfg["channels"]
    fs = p["eo_mi"]["fs"]
    data = io.load_edf_by_index(subject_file(cfg, "eo_mi", sid), chans["eo_mi_index"],
                                chans["canonical"], fs_expected=fs)
    data = _bandpass(cfg, band, data, fs)                    # whole recording
    data = segment.first_seconds(data, fs, p["eo_mi"]["duration_s"])
    base = {"dataset": "eo_mi", "subject": sid, "severity": "-", "lesion": "-"}
    return _rows(base, "EO-MI", _healthy_pieces(cfg, data, mode), chans["canonical"],
                 p["pe"]["mu"], p["pe"]["tau"])


def _ec_rs(cfg, sid, mode, band):
    p, chans = cfg["params"], cfg["channels"]
    fs = p["ec_rs"]["fs"]
    epochs = io.read_set_epochs(subject_file(cfg, "ec_rs", sid))
    epochs = io.select_by_index(epochs, chans["ec_rs_index"], chans["canonical"])
    epochs = _bandpass(cfg, band, epochs, fs, axis=1)        # each epoch separately (axis 1 = samples)
    data = segment.first_seconds(io.concat_epochs(epochs), fs, p["ec_rs"]["duration_s"])
    base = {"dataset": "ec_rs", "subject": sid, "severity": "-", "lesion": "-"}
    return _rows(base, "EC-RS", _healthy_pieces(cfg, data, mode), chans["canonical"],
                 p["pe"]["mu"], p["pe"]["tau"])


_WORKERS = {"stroke": _stroke, "eo_mi": _eo_mi, "ec_rs": _ec_rs}


def _task(args):
    cfg, dataset, sid, mode, band = args
    return _WORKERS[dataset](cfg, sid, mode, band)


def compute_pe_table(cfg, mode="paper", band="broadband", jobs=1, datasets=("stroke", "eo_mi", "ec_rs")):
    if mode not in ("paper", "published"):
        raise ValueError("mode must be 'paper' or 'published'")
    tasks = [(cfg, ds, sid, mode, band) for ds in datasets for sid in cfg["subjects"][ds]["ids"]]
    rows = []
    if jobs and jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for (_, ds, sid, _, _), r in zip(tasks, ex.map(_task, tasks)):
                print(f"  {ds:7s} {sid}  done", flush=True)
                rows += r
    else:
        for t in tasks:
            print(f"  {t[1]:7s} {t[2]}", flush=True)
            rows += _task(t)
    return pd.DataFrame(rows, columns=COLUMNS)
