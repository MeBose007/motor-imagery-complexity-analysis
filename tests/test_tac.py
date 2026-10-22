import numpy as np
from strokepe.tac import tac


def authors_tac(ts, mu, tau):
    pats = []
    for i in range(len(ts) - (mu - 1) * tau):
        v = np.array([ts[i + j * tau] for j in range(mu)], dtype=float)
        pats.append(tuple(np.argsort(np.argsort(v, kind="stable"), kind="stable")))
    uniq = sorted(set(pats))
    sym = np.array([uniq.index(p) for p in pats])
    n = sym.max() + 1
    m = np.zeros((n, n))
    for a, b in zip(sym[:-1], sym[1:]):
        m[a, b] += 1
    rs = m.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1
    m /= rs
    num = sum(abs(m[a, b] - m[b, a]) for a in range(n) for b in range(n) if a != b)
    den = sum(m[a, b] + m[b, a] for a in range(n) for b in range(n) if a != b)
    return 0.0 if den == 0 else num / den


def test_matches_authors_implementation():
    rng = np.random.default_rng(0)
    for mu in (2, 3):
        for _ in range(5):
            x = np.cumsum(rng.standard_normal(500))
            assert abs(tac(x, mu, 1) - authors_tac(x, mu, 1)) < 1e-12


def test_sine_is_fully_asymmetric():
    x = np.sin(2 * np.pi * np.arange(1300) / 1300)
    assert tac(x, 3, 1) > 0.99
