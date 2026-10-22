"""Transition asymmetry coefficient a_c (Supplementary Sec. I B, Fig. S2)."""
import numpy as np
from .pe import ordinal_codes


def transition_matrix(symbols):
    """Row-normalised transition matrix M_ab = P(next pattern = b | pattern = a),
    over the patterns that occur in ``symbols``."""
    _, s = np.unique(np.asarray(symbols), return_inverse=True)
    n = s.max() + 1
    counts = np.zeros((n, n))
    np.add.at(counts, (s[:-1], s[1:]), 1)
    rows = counts.sum(axis=1, keepdims=True)
    rows[rows == 0] = 1.0
    return counts / rows


def asymmetry_coefficient(symbols):
    """a_c = sum_{a != b} |M_ab - M_ba| / sum_{a != b} (M_ab + M_ba)."""
    m = transition_matrix(symbols)
    off = ~np.eye(m.shape[0], dtype=bool)
    num = np.abs(m - m.T)[off].sum()
    den = (m + m.T)[off].sum()
    return 0.0 if den == 0 else float(num / den)


def tac(x, mu=3, tau=1):
    """a_c of the ordinal-pattern sequence of ``x``."""
    return asymmetry_coefficient(ordinal_codes(x, mu, tau))


def windowed_tac(x, fs, window_s=1.0, mu=3, tau=1):
    """a_c in consecutive non-overlapping windows of ``window_s`` seconds."""
    w = int(round(window_s * fs))
    n = len(x) // w
    return np.array([tac(x[i * w:(i + 1) * w], mu, tau) for i in range(n)])
