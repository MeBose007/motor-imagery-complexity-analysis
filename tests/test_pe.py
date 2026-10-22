import math
from collections import Counter
import numpy as np
from strokepe.pe import permutation_entropy, ordinal_codes


def authors_pe(ts, mu, tau):
    """The authors' implementation (main analysis scripts)."""
    pats = []
    for i in range(len(ts) - (mu - 1) * tau):
        v = np.array([ts[i + j * tau] for j in range(mu)], dtype=float)
        pats.append(tuple(np.argsort(np.argsort(v, kind="stable"), kind="stable")))
    p = np.array(list(Counter(pats).values())) / len(pats)
    return -np.sum(p * np.log2(p)) / (math.lgamma(mu + 1) / np.log(2))


def test_matches_authors_implementation_including_ties():
    rng = np.random.default_rng(1)
    for mu in (3, 4, 5):
        for x in (rng.standard_normal(3000), rng.integers(0, 4, 3000).astype(float)):
            assert abs(permutation_entropy(x, mu, 1) - authors_pe(x, mu, 1)) < 1e-12
    x = rng.standard_normal(2000)
    assert abs(permutation_entropy(x, 5, 3) - authors_pe(x, 5, 3)) < 1e-12


def test_limits():
    assert permutation_entropy(np.arange(1000.0), 5, 1) == 0.0
    assert permutation_entropy(np.random.default_rng(2).standard_normal(200000), 5, 1) > 0.99


def test_codes_are_permutations():
    codes = ordinal_codes(np.random.default_rng(3).standard_normal(50000), 3, 1)
    assert len(np.unique(codes)) == 6
