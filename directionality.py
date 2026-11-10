"""
directionality.py - can GENIE3 recover the DIRECTION of true edges?

For every true edge i->j (whose reverse j->i is NOT a true edge), compare the
score W[i,j] with W[j,i]. Chance = 50% (ties excluded). Pearson |r| is symmetric,
so every pair is a tie: it cannot orient edges.

CAVEAT: the pooled data contain knockouts/knockdowns (interventions). Perturbing
gene i changes its targets but not its regulators, which by itself carries
directional information. So a good result here may come partly from the data
type, not only from the tree-based method.
"""
import numpy as np, pandas as pd
from scipy.stats import binomtest
from genie3 import genie3, pearson_baseline
from sweep import load_pooled


def direction_stats(W, gold):
    pos = set(zip(gold[gold.label == 1].regulator, gold[gold.label == 1].target))
    wins = losses = ties = 0
    for r, t in pos:
        if (t, r) in pos:          # reciprocal edge: direction is ambiguous, skip
            continue
        a, b = W.loc[r, t], W.loc[t, r]
        if a > b: wins += 1
        elif a < b: losses += 1
        else: ties += 1
    n = wins + losses
    p = binomtest(wins, n, 0.5).pvalue if n else float("nan")
    return wins, losses, ties, p


if __name__ == "__main__":
    rows = []
    for net in range(1, 6):
        expr, gold = load_pooled(net)
        for name, W in (("GENIE3 RF", genie3(expr, n_trees=200, max_features=0.1, method="RF")),
                        ("Pearson", pearson_baseline(expr))):
            w, l, t, p = direction_stats(W, gold)
            rows.append({"network": net, "method": name, "correct": w, "wrong": l, "ties": t,
                         "frac_correct_of_decided": round(w / (w + l), 3) if w + l else np.nan,
                         "binom_p": p})
            print(rows[-1], flush=True)
    df = pd.DataFrame(rows); df.to_csv("results/directionality.csv", index=False)
    print(df.to_string(index=False))
