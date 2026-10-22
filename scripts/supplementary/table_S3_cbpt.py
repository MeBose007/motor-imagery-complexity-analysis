"""Table S3 - cluster-based permutation test (Maris & Oostenveld) between every pair of
cohorts on per-segment PE (same pooling as Tables S1-S2), k-nearest-neighbour sensor
adjacency from the electrode positions of the stroke dataset, two-sided Student's t,
cluster mass = sum of t. Needs: 01_compute_pe.py --mode published, and electrodes_tsv in
config/paths.yaml. p-values have a floor of 1/(N+1)."""
import sys, pathlib; sys.path.insert(0, str(next(p for p in pathlib.Path(__file__).resolve().parents if (p / "strokepe").is_dir())))
import pandas as pd
from strokepe.cli import parser, setup, load_pe_table, markdown
from strokepe import stats, cbpt

COMPARISONS = [("minor", "severe", "M-S"), ("minor", "EO", "M-EO"), ("minor", "EC", "M-EC"),
               ("severe", "EO", "S-EO"), ("severe", "EC", "S-EC"), ("EO", "EC", "EO-EC")]
SEED_KEY = {"M-S": "minor-severe", "M-EO": "minor-EO", "S-EO": "severe-EO", "EO-EC": "EO-EC",
            "M-EC": "minor-EC", "S-EC": "severe-EC"}


def main(argv=None):
    p = parser(__doc__, mode=False)
    prm = None
    p.add_argument("--n-perm", type=int, nargs="+", default=None)
    p.add_argument("--k", type=int, nargs="+", default=None)
    p.add_argument("--alpha-form", type=float, default=None)
    args = p.parse_args(argv)
    cfg, out = setup(args)
    prm = cfg["params"]["cbpt"]
    n_perms, ks = args.n_perm or prm["n_perm"], args.k or prm["k_neighbors"]
    alpha = args.alpha_form or prm["alpha_form"]
    out = out / "supplementary"
    out.mkdir(exist_ok=True)

    frame = stats.trial_level_frame(load_pe_table(cfg, "published", args.band))
    channels = sorted(cfg["channels"]["canonical"])
    xyz = cbpt.load_electrodes(cfg["paths"]["electrodes_tsv"], channels)
    order = {"minor": cfg["subjects"]["stroke"]["ids"], "severe": cfg["subjects"]["stroke"]["ids"],
             "EO": cfg["subjects"]["eo_mi"]["ids"], "EC": cfg["subjects"]["ec_rs"]["ids"]}
    mats = {c: cbpt.cohort_matrix(frame, c, channels, order[c]) for c in order}
    rows = []
    for n_perm in n_perms:
        for k in ks:
            adj = cbpt.knn_adjacency(xyz, k)
            for a, b, lab in COMPARISONS:
                r = cbpt.cluster_permutation_test(mats[a], mats[b], adj, n_perm, alpha,
                                                  seed=prm["seeds"][SEED_KEY[lab]])
                sig = [c for c in r["clusters"] if c["p"] < 0.05]
                rows.append({"N": n_perm, "alpha": alpha, "k": k, "comparison": lab, "p": r["p"],
                             "significant clusters": " | ".join(
                                 f"{','.join(channels[i] for i in c['channels'])} (p={c['p']:.2g})" for c in sig)})
                print(f"N={n_perm} k={k} {lab:6s} p={r['p']:.3g}  {rows[-1]['significant clusters']}", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(out / "tableS3_cbpt.csv", index=False)
    (out / "tableS3_cbpt.md").write_text("# Table S3 - cluster-based permutation test\n\n"
                                         + markdown(res, floatfmt="{:.3g}") + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
