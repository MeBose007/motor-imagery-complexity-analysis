"""Ordinal patterns and permutation entropy (Bandt & Pompe 2002), Eqs. (1)-(5).

Ties are broken by temporal order (the earlier sample gets the lower rank), i.e. a
"predefined tie-breaking rule" as stated in Sec. II B. This is identical to the
authors' implementation ``np.argsort(np.argsort(x, kind="stable"), kind="stable")``.
"""
import math
import numpy as np


def ordinal_codes(x, mu=5, tau=1):
    """Integer code of the ordinal pattern of every embedding vector of ``x``."""
    x = np.asarray(x, dtype=float)
    n = x.size - (mu - 1) * tau
    if n <= 0:
        raise ValueError(f"series of length {x.size} too short for mu={mu}, tau={tau}")
    idx = np.arange(n)[:, None] + tau * np.arange(mu)[None, :]
    emb = x[idx]
    ranks = np.argsort(np.argsort(emb, axis=1, kind="stable"), axis=1, kind="stable")
    weights = mu ** np.arange(mu - 1, -1, -1)
    return ranks @ weights


def pattern_probabilities(x, mu=5, tau=1):
    """Relative frequency of each observed ordinal pattern, Eq. (3)."""
    _, counts = np.unique(ordinal_codes(x, mu, tau), return_counts=True)
    return counts / counts.sum()


def permutation_entropy(x, mu=5, tau=1, normalize=True):
    """Shannon entropy of the ordinal-pattern distribution, Eq. (4); divided by
    log(mu!) when ``normalize`` is True, Eq. (5)."""
    p = pattern_probabilities(x, mu, tau)
    h = float(-np.sum(p * np.log2(p)))
    if normalize:
        h /= math.lgamma(mu + 1) / math.log(2)
    return h
