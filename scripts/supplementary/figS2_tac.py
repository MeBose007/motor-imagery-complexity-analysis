"""Fig. S2 - time-resolved transition asymmetry coefficient a_c (mu = 2 and 3, tau = 1)
in non-overlapping 1-s windows over 80 s, averaged over channels and then over subjects,
with a periodic signal and pink noise as references."""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from strokepe.cli import parser, setup, save_figure
from strokepe.config import subject_file, severity_of
from strokepe import io, segment, synthetic
from strokepe.tac import windowed_tac

COLORS = {"Minor L-MIT": "#FF6B6B", "Minor R-MIT": "#FFA07A", "Severe L-MIT": "#8B0000",
          "Severe R-MIT": "#DC143C", "EO-MI": "#069AF3", "EC-RS": "#A9A9A9",
          "Periodic signal": "gold", "Pink noise": "#9B59B6", "Random noise": "green"}


def channel_mean_tac(data, fs, window_s, mu, tau):
    return np.mean([windowed_tac(ch, fs, window_s, mu, tau) for ch in data], axis=0)


def cohort_curves(cfg, mu, tau, window_s):
    p, chans = cfg["params"], cfg["channels"]
    canonical = chans["canonical"]
    curves = {k: [] for k in ["Minor L-MIT", "Minor R-MIT", "Severe L-MIT", "Severe R-MIT", "EO-MI", "EC-RS"]}
    for sid in cfg["subjects"]["stroke"]["ids"]:
        data = io.load_stroke_csv(subject_file(cfg, "stroke", sid), chans["stroke_drop_columns"], canonical)
        segs = segment.stroke_mi_segments(data, p["stroke"]["trial_samples"], p["stroke"]["mi_start"],
                                          p["stroke"]["mi_stop"], p["stroke"]["trials_expected"])
        sev = severity_of(cfg, sid).capitalize()
        for cond in ("L-MIT", "R-MIT"):
            x = np.concatenate([s for _, c, s in segs if c == cond], axis=1)
            curves[f"{sev} {cond}"].append(channel_mean_tac(x, p["stroke"]["fs"], window_s, mu, tau))
        print(f"  stroke {sid}", flush=True)
    for sid in cfg["subjects"]["eo_mi"]["ids"]:
        fs = p["eo_mi"]["fs"]
        x = io.load_edf_by_index(subject_file(cfg, "eo_mi", sid), chans["eo_mi_index"], canonical, fs)
        curves["EO-MI"].append(channel_mean_tac(segment.first_seconds(x, fs, p["eo_mi"]["duration_s"]),
                                                fs, window_s, mu, tau))
    for sid in cfg["subjects"]["ec_rs"]["ids"]:
        fs = p["ec_rs"]["fs"]
        ep = io.select_by_index(io.read_set_epochs(subject_file(cfg, "ec_rs", sid)), chans["ec_rs_index"], canonical)
        x = segment.first_seconds(io.concat_epochs(ep), fs, p["ec_rs"]["duration_s"])
        curves["EC-RS"].append(channel_mean_tac(x, fs, window_s, mu, tau))
    return {k: np.mean(v, axis=0) for k, v in curves.items() if v}


def reference_curves(cfg, mu, tau, window_s, white):
    t = cfg["params"]["tac"]
    fs, dur = t["synthetic_fs"], 80
    rng = np.random.default_rng(t["seed"])
    sigs = {"Periodic signal": synthetic.sine(fs, dur, t["sine_hz"], t["sine_amplitude"]),
            "Pink noise": synthetic.pink_noise(fs, dur, rng)}
    if white:
        sigs["Random noise"] = synthetic.white_noise(fs, dur, rng)
    return {k: windowed_tac(v, fs, window_s, mu, tau) for k, v in sigs.items()}


def main(argv=None):
    p = parser(__doc__, mode=False, band=False)
    p.add_argument("--mu", type=int, nargs="+", default=[2, 3])
    p.add_argument("--include-white-noise", action="store_true",
                   help="also plot white noise (not shown in the published figure)")
    args = p.parse_args(argv)
    cfg, _ = setup(args)
    out = setup(args)[1].parent / "supplementary"
    out.mkdir(parents=True, exist_ok=True)
    tau, window_s = cfg["params"]["pe"]["tau"], cfg["params"]["tac"]["window_s"]
    fig, axes = plt.subplots(1, len(args.mu), figsize=(6 * len(args.mu), 6), sharey=True, squeeze=False)
    for ax, mu, lab in zip(axes[0], args.mu, "abcdef"):
        print(f"mu = {mu}")
        curves = {**cohort_curves(cfg, mu, tau, window_s),
                  **reference_curves(cfg, mu, tau, window_s, args.include_white_noise)}
        pd.DataFrame({k: pd.Series(v) for k, v in curves.items()}).rename_axis("window").to_csv(out / f"figureS2_tac_mu{mu}.csv", float_format="%.6f")
        for name, y in curves.items():
            ls = "--" if name in ("Periodic signal", "Pink noise", "Random noise") else "-"
            ax.plot(np.arange(len(y)) * window_s, y, ls, color=COLORS[name], lw=2, label=name)
        ax.set_title(f"({lab})", loc="left", fontsize=18)
        ax.text(0.03, 0.95, f"$\\mu$={mu}", transform=ax.transAxes, va="top", fontsize=14)
        ax.set_xlabel("Time (s)", fontsize=18)
        ax.set_ylim(0, 1.01)
    axes[0][0].set_ylabel("$a_c$", fontsize=20)
    axes[0][0].legend(fontsize=9, loc="center right")
    fig.tight_layout()
    save_figure(fig, out, "figureS2")


if __name__ == "__main__":
    main()
