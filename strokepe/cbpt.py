"""Cluster-based permutation test (Maris & Oostenveld 2007), Table S3."""
from collections import deque
import numpy as np
import pandas as pd
from scipy import stats


def load_electrodes(tsv_path, channels):
    """3-D positions of ``channels`` from a BIDS electrodes.tsv (columns name|channel, x, y, z)."""
    t = pd.read_csv(tsv_path, sep="\t")
    cols = {c.lower(): c for c in t.columns}
    name_col = cols.get("channel", cols.get("name"))
    if name_col is None or not {"x", "y", "z"} <= set(cols):
        raise ValueError(f"{tsv_path}: need columns name (or channel), x, y, z; found {list(t.columns)}")
    key = lambda s: "".join(ch for ch in str(s) if ch.isalnum()).lower()
    pos = {key(r[name_col]): (r[cols["x"]], r[cols["y"]], r[cols["z"]]) for _, r in t.iterrows()}
    missing = [c for c in channels if key(c) not in pos]
    if missing:
        raise ValueError(f"{tsv_path}: no coordinates for {missing}")
    return np.array([pos[key(c)] for c in channels], dtype=float)


def knn_adjacency(xyz, k=2):
    """Symmetric k-nearest-neighbour adjacency from 3-D Euclidean distances."""
    d = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(axis=2))
    n = len(xyz)
    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        adj[i, np.argsort(d[i])[1:k + 1]] = True
    adj |= adj.T
    np.fill_diagonal(adj, False)
    return adj


def find_clusters(mask, adj):
    seen = np.zeros(len(mask), dtype=bool)
    clusters = []
    for i in range(len(mask)):
        if mask[i] and not seen[i]:
            seen[i] = True
            q, cl = deque([i]), [i]
            while q:
                u = q.popleft()
                for v in np.where(adj[u])[0]:
                    if mask[v] and not seen[v]:
                        seen[v] = True
                        q.append(v)
                        cl.append(v)
            clusters.append(cl)
    return clusters


def _t_student(xa, xb):
    """Student's t (equal variances) along axis -2 for stacked samples (..., n, channels)."""
    na, nb = xa.shape[-2], xb.shape[-2]
    ma, mb = xa.mean(axis=-2), xb.mean(axis=-2)
    va, vb = xa.var(axis=-2, ddof=1), xb.var(axis=-2, ddof=1)
    sp = ((na - 1) * va + (nb - 1) * vb) / (na + nb - 2)
    return (ma - mb) / np.sqrt(sp * (1 / na + 1 / nb))


def _max_mass(t, tcrit, adj):
    cl = find_clusters(np.abs(t) > tcrit, adj)
    return 0.0 if not cl else max(abs(t[c].sum()) for c in cl)


def cluster_permutation_test(xa, xb, adj, n_perm=10000, alpha_form=0.05, seed=0, batch=200):
    """Two-sided test of xa vs xb (observations x channels). Cluster statistic = sum of t
    within a cluster (cluster mass); test statistic = largest |mass|. The permutation
    sequence (one rng.permutation per iteration) matches the authors' implementation."""
    x = np.vstack([xa, xb])
    n0, ntot = len(xa), len(xa) + len(xb)
    t_obs = stats.ttest_ind(xa, xb, axis=0, equal_var=True).statistic
    tcrit = stats.t.ppf(1 - alpha_form / 2, ntot - 2)
    clusters = find_clusters(np.abs(t_obs) > tcrit, adj)
    masses = [float(t_obs[c].sum()) for c in clusters]
    s_obs = max((abs(m) for m in masses), default=0.0)

    rng = np.random.default_rng(seed)
    perm_max = np.empty(n_perm)
    done = 0
    while done < n_perm:
        b = min(batch, n_perm - done)
        perms = np.stack([rng.permutation(ntot) for _ in range(b)])
        xp = x[perms]                                   # (b, ntot, channels)
        t_perm = _t_student(xp[:, :n0], xp[:, n0:])
        perm_max[done:done + b] = [_max_mass(t, tcrit, adj) for t in t_perm]
        done += b

    p_max = (np.sum(perm_max >= s_obs) + 1) / (n_perm + 1)
    reports = sorted(
        [{"channels": c, "mass": m, "p": (np.sum(perm_max >= abs(m)) + 1) / (n_perm + 1)}
         for c, m in zip(clusters, masses)], key=lambda r: abs(r["mass"]), reverse=True)
    return {"p": p_max, "tcrit": tcrit, "t_obs": t_obs, "clusters": reports, "perm_max": perm_max}


def cohort_matrix(frame, cohort, channels, subject_order):
    """Observations x channels for one cohort; rows ordered trial-first, then subject,
    as in the authors' concatenation."""
    s = frame[frame["cohort"] == cohort]
    wide = s.pivot_table(index=["trial", "subject"], columns="channel", values="pe")
    order = {sid: i for i, sid in enumerate(subject_order)}
    idx = sorted(wide.index, key=lambda ts: (ts[0], order.get(ts[1], 10**9)))
    return wide.loc[idx, channels].to_numpy()
