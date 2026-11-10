"""
evaluate.py - score an inferred edge-score matrix against the DREAM4 gold standard.
"""
import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score


def evaluate(W, gold):
    """W: DataFrame of edge scores W.loc[i, j] for i -> j.
    gold: DataFrame (regulator, target, label).
    Evaluates all ordered pairs (i != j). Returns dict with AUPR and AUROC."""
    scores, labels = [], []
    for r, t, l in zip(gold.regulator, gold.target, gold.label):
        if r == t:
            continue
        scores.append(W.loc[r, t])
        labels.append(l)
    scores, labels = np.array(scores), np.array(labels)
    return {"AUPR": average_precision_score(labels, scores),
            "AUROC": roc_auc_score(labels, scores),
            "n_pairs": len(labels), "n_true": int(labels.sum()),
            "random_AUPR": labels.mean()}   # AUPR expected by chance = edge density
