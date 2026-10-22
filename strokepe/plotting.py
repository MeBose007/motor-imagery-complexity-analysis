"""Figure helpers. Colours and styling follow the published figures."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

COLORS = {"L-MIT": "#8C000F", "R-MIT": "#7FFF00", "EO-MI": "#069AF3", "EC-RS": "lightgrey"}
EPS = r"$\mathcal{E}$"


def violin_box(ax, groups, colors, ylim=None, label_size=20):
    """Violin (KDE) with a narrow box plot inside; groups = {label: values}."""
    labels = list(groups)
    df = pd.DataFrame({"x": np.repeat(range(len(labels)), [len(groups[k]) for k in labels]),
                       "y": np.concatenate([np.asarray(groups[k], float) for k in labels])})
    pal = {i: colors[k] for i, k in enumerate(labels)}
    sns.violinplot(data=df, x="x", y="y", hue="x", palette=pal, inner=None, linewidth=1.5,
                   alpha=0.7, legend=False, ax=ax)
    for i, k in enumerate(labels):
        ax.boxplot(np.asarray(groups[k], float), positions=[i], widths=0.1, patch_artist=True,
                   boxprops=dict(facecolor="white", edgecolor="black"),
                   whiskerprops=dict(color="black", linewidth=1.5),
                   capprops=dict(color="black", linewidth=1.5),
                   medianprops=dict(color="red", linewidth=1.5),
                   flierprops=dict(marker="o", markerfacecolor="none",
                                   markeredgecolor=colors[k], markersize=4))
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([k.split("|")[-1] for k in labels])
    ax.set_xlabel("")
    ax.set_ylabel(EPS, fontsize=30)
    ax.tick_params(axis="both", labelsize=label_size)
    if ylim:
        ax.set_ylim(*ylim)


def significance_bars(ax, pairs, groups, step=0.035):
    """Draw horizontal bars with '*' (p < 0.05) or 'ns' between violins.
    pairs = [(i, j, p), ...] with x positions i, j."""
    top = max(np.max(v) for v in groups.values())
    span = ax.get_ylim()[1] - ax.get_ylim()[0]
    for n, (i, j, p) in enumerate(sorted(pairs, key=lambda q: abs(q[1] - q[0]))):
        y = top + span * (0.04 + step * n)
        ax.plot([i, i, j, j], [y - span * 0.01, y, y, y - span * 0.01], color="black", lw=1.2)
        ax.text((i + j) / 2, y, "*" if p < 0.05 else "ns", ha="center",
                va="bottom" if p >= 0.05 else "center", fontsize=16 if p < 0.05 else 11)


def group_brackets(ax, spans, fontsize=12):
    """Brackets with a label under pairs of violins: spans = [(start, end, text), ...]."""
    y0, y1 = ax.get_ylim()
    r = y1 - y0
    for a, b, text in spans:
        y = y0 - 0.20 * r
        ax.plot([a, a, b, b], [y + 0.02 * r, y, y, y + 0.02 * r], "k-", lw=2, clip_on=False)
        ax.text((a + b) / 2, y - 0.03 * r, text, ha="center", va="top", fontsize=fontsize, clip_on=False)


def channel_boxplots(ax, data, channels, conditions, ylim=(0.2, 0.8)):
    """data[(condition, channel)] -> values. Four boxes per channel (Fig. 3a,b)."""
    width = len(conditions) + 1
    for ci, ch in enumerate(channels):
        for k, cond in enumerate(conditions):
            vals = np.asarray(data[(cond, ch)], float)
            ax.boxplot(vals, positions=[ci * width + k], widths=0.7, patch_artist=True,
                       boxprops=dict(facecolor=COLORS[cond], alpha=0.8),
                       medianprops=dict(color="black", linewidth=1.5),
                       flierprops=dict(marker="o", markersize=3, markerfacecolor="none",
                                       markeredgecolor=COLORS[cond]))
        if ci:
            ax.axvline(ci * width - 1, color="grey", ls=":", lw=0.8)
    ax.set_xticks([ci * width + (len(conditions) - 1) / 2 for ci in range(len(channels))])
    ax.set_xticklabels(channels, fontsize=13)
    ax.set_xlim(-1, len(channels) * width - 1)
    ax.set_ylim(*ylim)
    ax.set_ylabel(EPS, fontsize=25)
    ax.tick_params(axis="y", labelsize=14)
    handles = [plt.Rectangle((0, 0), 1, 1, color=COLORS[c]) for c in conditions]
    ax.legend(handles, conditions, loc="upper right", fontsize=9, frameon=False)


# ------------------------------------------------------------------ topographic maps
_MANUAL_POS = {  # used only if a channel is missing from the MNE standard montages
    "FCz": [0.0, 0.094, 0.044], "FC3": [-0.055, 0.067, 0.051], "FC4": [0.055, 0.067, 0.051],
    "FT7": [-0.082, 0.026, -0.009], "FT8": [0.082, 0.026, -0.009],
    "CP3": [-0.055, -0.067, 0.051], "CP4": [0.055, -0.067, 0.051],
    "TP7": [-0.082, -0.026, -0.009], "TP8": [0.082, -0.026, -0.009]}


def montage_info(ch_names):
    """MNE Info with the authors' montage: standard_1005/1020 positions shifted up by 1 cm."""
    import warnings
    import mne
    pos = {}
    for name, new_name in (("standard_1005", "colin27_1005"), ("standard_1020", "colin27_1020")):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)   # renamed in MNE 1.14
            try:
                m = mne.channels.make_standard_montage(name)
            except ValueError:
                m = mne.channels.make_standard_montage(new_name)
        pos.update(m.get_positions()["ch_pos"])
    for ch, p in _MANUAL_POS.items():
        pos.setdefault(ch, np.array(p))
    ch_pos = {}
    for ch in ch_names:
        p = np.array(pos.get(ch, [0.0, 0.0, 0.1]), float).copy()
        p[1] += 0.01
        ch_pos[ch] = p
    info = mne.create_info(ch_names=list(ch_names), sfreq=100, ch_types="eeg")
    info.set_montage(mne.channels.make_dig_montage(ch_pos=ch_pos, coord_frame="head"))
    return info


def topoplot(ax, values, info):
    """Min-max normalised map (each map scaled to [0, 1] separately), as in Figs. 3-4."""
    import mne
    v = np.asarray(values, float)
    norm = (v - v.min()) / (v.max() - v.min())
    im, _ = mne.viz.plot_topomap(norm, info, axes=ax, show=False, sphere=(0., 0., 0., 0.18),
                                 contours=8, cmap="RdBu_r", sensors=False, res=64,
                                 outlines="head", vlim=(0, 1))
    ax.axis("off")
    return im
