"""Write a small synthetic dataset in the same file formats as the real data, together
with a reduced configuration, so the whole pipeline can be exercised without the EEG data.
The signals are random and carry no meaning."""
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
STROKE_IDS = {"minor": ["sub-02", "sub-03", "sub-07", "sub-20"], "severe": ["sub-01", "sub-09", "sub-11", "sub-15"]}
LESION = {"LHA": ["sub-01", "sub-07", "sub-09", "sub-20"], "RHA": ["sub-02", "sub-03", "sub-11", "sub-15"]}
EO_IDS = ["S017", "S022", "S028", "S030"]
EC_IDS = ["sub-003", "sub-020", "sub-021", "sub-026"]
DORTMUND_IDS = ["sub-002", "sub-014"]
EXTRA_STROKE_COLS = {15: "T3", 16: "T4", 17: "CPz", 25: "T5", 26: "T6"}


def _noise(rng, n_ch, n, smooth):
    x = rng.standard_normal((n_ch, n + smooth))
    k = np.ones(smooth) / smooth
    return np.array([np.convolve(r, k, mode="valid")[:n] for r in x])


def _write_edf(path, data, fs, labels):
    import pyedflib
    f = pyedflib.EdfWriter(str(path), len(labels), file_type=pyedflib.FILETYPE_EDFPLUS)
    try:
        hdrs = []
        for lab, ch in zip(labels, data):
            lim = float(np.ceil(np.abs(ch).max() + 1))
            hdrs.append({"label": lab, "dimension": "uV", "sample_frequency": fs,
                         "physical_min": -lim, "physical_max": lim,
                         "digital_min": -32768, "digital_max": 32767})
        f.setSignalHeaders(hdrs)
        f.writeSamples(list(data))
    finally:
        f.close()


def make(tmp):
    tmp = Path(tmp)
    rng = np.random.default_rng(0)
    cfgdir, data = tmp / "config", tmp / "data"
    for d in ["stroke", "eo_mi", "ec_rs", "dortmund"]:
        (data / d).mkdir(parents=True, exist_ok=True)
    cfgdir.mkdir(exist_ok=True)
    channels = yaml.safe_load(open(ROOT / "config/channels.yaml"))
    params = yaml.safe_load(open(ROOT / "config/params.yaml"))
    subjects = yaml.safe_load(open(ROOT / "config/subjects.yaml"))
    canon = channels["canonical"]

    # reduced sizes keep the test fast
    params["stroke"].update(trial_samples=800, mi_start=200, mi_stop=600)
    params["eo_mi"]["duration_s"] = 20
    params["ec_rs"]["duration_s"] = 20
    params["psd"]["duration_s"] = 20
    ids = STROKE_IDS["minor"] + STROKE_IDS["severe"]
    subjects["stroke"].update(ids=sorted(ids), severity=STROKE_IDS, lesion_side=LESION)
    subjects["eo_mi"]["ids"] = EO_IDS
    subjects["ec_rs"]["ids"] = EC_IDS
    subjects["dortmund"]["ids"] = DORTMUND_IDS
    for name, obj in [("channels", channels), ("params", params), ("subjects", subjects)]:
        yaml.safe_dump(obj, open(cfgdir / f"{name}.yaml", "w"), sort_keys=False)
    yaml.safe_dump({"stroke_dir": str(data / "stroke"), "eo_mi_dir": str(data / "eo_mi"),
                    "ec_rs_dir": str(data / "ec_rs"), "dortmund_dir": str(data / "dortmund"),
                    "electrodes_tsv": str(data / "stroke" / "task-motor-imagery_electrodes.tsv"),
                    "results_dir": str(tmp / "results")}, open(cfgdir / "paths.yaml", "w"))

    # stroke CSVs: 30 columns, 5 extra channels at the positions removed by the pipeline
    names, it = [], iter(canon)
    for i in range(30):
        names.append(EXTRA_STROKE_COLS.get(i) or next(it))
    n = 40 * params["stroke"]["trial_samples"]
    for sid in subjects["stroke"]["ids"]:
        x = _noise(rng, 30, n, smooth=4 if sid in STROKE_IDS["minor"] else 5) * 10
        pd.DataFrame(x.T, columns=names).to_csv(
            data / "stroke" / subjects["stroke"]["file_pattern"].format(id=sid), index=False, float_format="%.5f")

    # PhysioNet-like EDF: 64 channels at 160 Hz (labels arbitrary, channels picked by index)
    for sid in EO_IDS:
        x = _noise(rng, 64, 160 * 25, smooth=1) * 20
        _write_edf(data / "eo_mi" / subjects["eo_mi"]["file_pattern"].format(id=sid), x, 160,
                   [f"E{i:02d}" for i in range(64)])

    # SRM-like .set: 64 channels, 1024 Hz, 2-s epochs
    import scipy.io
    for sid in EC_IDS:
        t = np.arange(1024 * 2 * 12) / 1024
        x = _noise(rng, 64, t.size, smooth=12) * 5 + 20 * np.sin(2 * np.pi * 10 * t)
        epochs = x.reshape(64, 12, 2048).transpose(0, 2, 1)
        scipy.io.savemat(data / "ec_rs" / subjects["ec_rs"]["file_pattern"].format(id=sid), {"data": epochs})

    # Dortmund-like EDF with channel names, 1000 Hz
    labels = channels["psd_left"] + channels["psd_right"] + ["Cz", "Pz"]
    for sid in DORTMUND_IDS:
        for task in ("EyesClosed", "EyesOpen"):
            x = _noise(rng, len(labels), 1000 * 22, smooth=6) * 10
            _write_edf(data / "dortmund" / subjects["dortmund"]["file_pattern"].format(id=sid, task=task),
                       x, 1000, labels)

    # BIDS-style electrode positions
    from strokepe.plotting import montage_info
    info = montage_info(canon + ["Iz"])
    pos = {ch["ch_name"]: ch["loc"][:3] for ch in info["chs"]}
    pd.DataFrame([{"name": c, "x": pos[c][0], "y": pos[c][1], "z": pos[c][2]} for c in canon + ["Iz"]]).to_csv(
        data / "stroke" / "task-motor-imagery_electrodes.tsv", sep="\t", index=False)
    return cfgdir


if __name__ == "__main__":
    import sys
    d = make(Path(sys.argv[1] if len(sys.argv) > 1 else "synthetic"))
    print(f"synthetic data written; set STROKEPE_CONFIG_DIR={d}")
