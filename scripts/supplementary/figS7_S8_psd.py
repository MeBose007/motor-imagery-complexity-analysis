"""Figs. S7-S8 - power spectral density (Welch, nperseg = 1024) of the first 180 s of the
Dortmund resting-state recordings, eyes closed (S7) and eyes open (S8), all subjects
overlaid for 20 channels."""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
from strokepe.cli import parser, setup, save_figure
from strokepe.config import subject_file
from strokepe import io, segment

BANDS = [(1, 4, "purple", "\u03b4"), (4, 8, "cyan", "\u03b8"), (8, 12, "green", "\u03b1"),
         (13, 30, "orange", "\u03b2"), (30, 50, "red", "\u03b3")]


def figure(cfg, task, out, name):
    chans = cfg["channels"]
    names = chans["psd_left"] + chans["psd_right"]
    prm = cfg["params"]["psd"]
    fig, axes = plt.subplots(len(chans["psd_left"]), 2, figsize=(14, 22), sharex=True)
    ids = cfg["subjects"]["dortmund"]["ids"]
    cmap = plt.get_cmap("tab20")
    for si, sid in enumerate(ids):
        data, fs = io.load_edf_by_name(subject_file(cfg, "dortmund", sid, task=task), names)
        data = segment.first_seconds(data, fs, prm["duration_s"])
        for ci, ch in enumerate(names):
            ax = axes[ci % len(chans["psd_left"]), ci // len(chans["psd_left"])]
            f, pxx = welch(data[ci], fs=fs, nperseg=min(prm["nperseg"], data.shape[1] // 4))
            ax.plot(f, 10 * np.log10(pxx), lw=1, alpha=0.7, color=cmap(si % 20), label=f"S{si + 1:02d}")
    for ci, ch in enumerate(names):
        ax = axes[ci % len(chans["psd_left"]), ci // len(chans["psd_left"])]
        for lo, hi, col, sym in BANDS:
            ax.axvspan(lo, hi, alpha=0.1, color=col)
        ax.set_xlim(0, 50)
        ax.text(0.98, 0.9, ch, transform=ax.transAxes, ha="right", fontsize=10,
                bbox=dict(facecolor="white", edgecolor="black", lw=0.5))
    fig.supxlabel("frequency (Hz)", fontsize=16)
    fig.supylabel("power (dB)", fontsize=16)
    axes[-1, 0].legend(ncol=9, fontsize=7, loc="upper center", bbox_to_anchor=(1.05, -0.35))
    fig.tight_layout()
    save_figure(fig, out, name)


def main(argv=None):
    p = parser(__doc__, mode=False, band=False)
    p.add_argument("--task", choices=["EyesClosed", "EyesOpen", "both"], default="both")
    args = p.parse_args(argv)
    cfg, out = setup(args)
    out = out.parent / "supplementary"
    out.mkdir(parents=True, exist_ok=True)
    if args.task in ("EyesClosed", "both"):
        figure(cfg, "EyesClosed", out, "figureS7_psd_eyes_closed")
    if args.task in ("EyesOpen", "both"):
        figure(cfg, "EyesOpen", out, "figureS8_psd_eyes_open")


if __name__ == "__main__":
    main()
