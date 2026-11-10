"""
sweep.py - experiments beyond the tutorial (the project's "change").

Data  : DREAM4 In Silico Size-100, pooled steady-state profiles
        (wild-type + single knockouts + single knockdowns = 201 samples x 100 genes).
        NOTE: the paper uses the *multifactorial* data. In the Bioconductor copy we
        downloaded, that file did not recover the gold standard (see README), so we
        use the verified steady-state files instead.

Experiments (all scored against the gold standard with AUPR / AUROC):
  A. Sweep of K (max_features) for Random Forest vs Extra-Trees
  B. Sweep of number of trees (Random Forest, default K)
  C. Seed repeats of the default setting (how much is noise?)
  D. Baseline: |Pearson correlation| (symmetric co-expression, from this course)
  --robust : default setting on all 5 networks (is network 4 typical?)

Usage:  python sweep.py --net 4            (full)
        python sweep.py --net 4 --quick    (fast smoke test)
        python sweep.py --robust
"""
import argparse, os, time
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from genie3 import genie3, pearson_baseline, load_expression, load_gold_standard
from evaluate import evaluate


def load_pooled(net, data_dir="data"):
    """Pool wild-type, knockout and knockdown steady states into one matrix."""
    d = f"{data_dir}/insilico_size100_{net}"
    parts = [load_expression(f"{d}/{f}.tsv") for f in ("wildtype", "knockouts", "knockdowns")]
    return pd.concat(parts, ignore_index=True), load_gold_standard(f"{d}/goldStandard.tsv")


def run(expr, gold, label, **kw):
    t0 = time.time()
    res = evaluate(genie3(expr, **kw), gold)
    out = {"experiment": label, "AUPR": res["AUPR"], "AUROC": res["AUROC"],
           "chance_AUPR": res["random_AUPR"], "seconds": round(time.time() - t0)}
    print(out, flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", type=int, default=4)
    ap.add_argument("--quick", action="store_true", help="tiny settings for a smoke test")
    ap.add_argument("--robust", action="store_true", help="default setting on all 5 networks")
    ap.add_argument("--data_dir", default="data")
    a = ap.parse_args()
    os.makedirs("results", exist_ok=True)

    # ---------- robustness across all five networks ----------
    if a.robust:
        rows = []
        n_trees = 20 if a.quick else 200
        for net in range(1, 6):
            expr, gold = load_pooled(net, a.data_dir)
            r = run(expr, gold, f"GENIE3 RF net{net}", n_trees=n_trees, max_features=0.1, method="RF")
            b = evaluate(pearson_baseline(expr), gold)
            rows.append({**r, "network": net, "method": "GENIE3",})
            rows.append({"experiment": f"Pearson net{net}", "AUPR": b["AUPR"], "AUROC": b["AUROC"],
                         "chance_AUPR": b["random_AUPR"], "seconds": 0, "network": net, "method": "Pearson"})
        df = pd.DataFrame(rows); df.to_csv("results/robust_all_networks.csv", index=False)
        print(df[["network", "method", "AUPR", "AUROC", "chance_AUPR"]].round(3).to_string(index=False))
        return

    expr, gold = load_pooled(a.net, a.data_dir)
    p = expr.shape[1]
    print(f"Network {a.net}: {expr.shape[0]} samples x {p} genes")
    T_K = 10 if a.quick else 100          # trees used in the K sweep
    rows = []

    # ---------- D. baseline ----------
    b = evaluate(pearson_baseline(expr), gold)
    base = {"experiment": "Pearson |r| baseline", "AUPR": b["AUPR"], "AUROC": b["AUROC"],
            "chance_AUPR": b["random_AUPR"], "seconds": 0}
    print(base); rows.append(base)

    # ---------- A. sweep K for RF and ET ----------
    k_fracs = [0.1, 0.33, 1.0]            # K ~ sqrt(p) (~10), p/3, p (all genes)
    for method in ("RF", "ET"):
        for f in k_fracs:
            K = max(1, int(round(f * (p - 1))))
            r = run(expr, gold, f"{method} K={K} trees={T_K}", n_trees=T_K, max_features=f, method=method)
            r.update(sweep="K", method=method, K=K, n_trees=T_K); rows.append(r)

    # ---------- B. sweep number of trees (RF, K=0.1) ----------
    for nt in ([5, 10, 20] if a.quick else [25, 100, 400]):
        r = run(expr, gold, f"RF trees={nt} K~sqrt(p)", n_trees=nt, max_features=0.1, method="RF")
        r.update(sweep="trees", method="RF", K=int(round(0.1 * (p - 1))), n_trees=nt); rows.append(r)

    # ---------- C. seed repeats of the default ----------
    seeds = [0, 1] if a.quick else [0, 1, 2]
    seed_rows = []
    for s in seeds:
        r = run(expr, gold, f"RF default seed={s}", n_trees=T_K, max_features=0.1, method="RF", random_state=s)
        r.update(sweep="seed", seed=s); rows.append(r); seed_rows.append(r)
    sa = np.array([r["AUPR"] for r in seed_rows])
    print(f"\nSeed repeat AUPR: mean={sa.mean():.4f} sd={sa.std(ddof=1) if len(sa)>1 else float('nan'):.4f}")

    df = pd.DataFrame(rows)
    df.to_csv(f"results/sweep_net{a.net}.csv", index=False)

    # ---------- plots ----------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    kd = df[df.sweep == "K"] if "sweep" in df else df.iloc[0:0]
    for m, c in (("RF", "tab:blue"), ("ET", "tab:orange")):
        s = kd[kd.method == m]; ax[0].plot(s.K, s.AUPR, "o-", label=m, color=c)
    ax[0].axhline(b["AUPR"], ls="--", color="gray", label="Pearson baseline")
    ax[0].axhline(b["random_AUPR"], ls=":", color="red", label="chance")
    ax[0].set_xlabel("K (candidate genes per split)"); ax[0].set_ylabel("AUPR")
    ax[0].set_title(f"Effect of K ({T_K} trees)"); ax[0].legend()
    td = df[df.sweep == "trees"] if "sweep" in df else df.iloc[0:0]
    ax[1].plot(td.n_trees, td.AUPR, "o-", color="tab:blue", label="GENIE3 RF")
    ax[1].axhline(b["AUPR"], ls="--", color="gray", label="Pearson baseline")
    ax[1].set_xscale("log"); ax[1].set_xlabel("number of trees"); ax[1].set_ylabel("AUPR")
    ax[1].set_title("Effect of number of trees"); ax[1].legend()
    plt.tight_layout(); plt.savefig(f"results/sweep_net{a.net}.png", dpi=150)
    print("\nSaved results/ (csv + png)")
    print(df[["experiment", "AUPR", "AUROC", "seconds"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
