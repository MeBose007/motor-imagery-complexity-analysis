"""Cutting recordings into the segments analysed in the paper."""
import warnings
import numpy as np


def stroke_mi_segments(data, trial_samples=4000, mi_start=1000, mi_stop=3000, trials_expected=40):
    """Split a continuous stroke recording (channels x samples) into trials and return
    the 4-s motor-imagery part of each one as (trial_index, condition, segment).
    Trials alternate left/right starting with left (trial 1 = L-MIT)."""
    n_trials = data.shape[1] // trial_samples
    if n_trials != trials_expected:
        warnings.warn(f"found {n_trials} trials, expected {trials_expected}")
    out = []
    for t in range(n_trials):
        start = t * trial_samples
        seg = data[:, start + mi_start:start + mi_stop]
        out.append((t, "L-MIT" if t % 2 == 0 else "R-MIT", seg))
    return out


def first_seconds(data, fs, seconds):
    n = int(round(fs * seconds))
    if data.shape[1] < n:
        warnings.warn(f"recording has {data.shape[1]} samples, fewer than {n} ({seconds} s)")
    return data[:, :n]


def equal_blocks(data, n_blocks=20):
    """Split channels x samples into ``n_blocks`` equal consecutive blocks
    (remainder samples at the end are dropped)."""
    size = data.shape[1] // n_blocks
    return [data[:, i * size:(i + 1) * size] for i in range(n_blocks)]
