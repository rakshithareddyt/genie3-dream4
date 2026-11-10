"""
genie3.py - GENIE3 (Huynh-Thu et al., 2010) re-implemented with scikit-learn.

Idea: for each target gene j, train a tree ensemble that predicts the expression
of gene j from ALL OTHER genes. The importance of gene i in that model becomes
the score of the directed edge  i -> j.  Stacking the p models gives a p x p
matrix of edge scores, which we rank to get a directed network.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor


def load_expression(path):
    """Read a DREAM4 TSV: rows = samples (conditions), columns = genes."""
    return pd.read_csv(path, sep="\t")


def load_gold_standard(path):
    """Read the DREAM4 gold standard: 3 columns (regulator, target, 0/1 label)."""
    return pd.read_csv(path, sep="\t", header=None,
                       names=["regulator", "target", "label"])


def _raw_importance(model, n_features):
    """Average UN-normalised impurity decrease over all trees.

    sklearn's .feature_importances_ rescales to sum to 1, which would hide how
    much variance the predictors actually explain. GENIE3 needs the raw
    variance reduction (our targets have variance 1, so it is comparable
    across genes), hence we compute it directly from each tree.
    """
    imp = np.zeros(n_features)
    for tree in model.estimators_:
        imp += tree.tree_.compute_feature_importances(normalize=False)
    return imp / len(model.estimators_)


def genie3(expr, n_trees=1000, max_features="sqrt", method="RF",
           random_state=0, n_jobs=-1):
    """Run GENIE3.

    expr         : DataFrame (samples x genes)
    n_trees      : trees per ensemble (paper uses 1000)
    max_features : K, number of candidate predictors tried at each split.
                   "sqrt" -> sqrt(p), 1.0 -> all p.
    method       : "RF" (Random Forest) or "ET" (Extra-Trees)
    Returns a DataFrame W where W.loc[i, j] = score of edge  i -> j.
    """
    genes = list(expr.columns)
    p = len(genes)
    # Standardise each gene to unit variance so importances are comparable.
    X_all = ((expr - expr.mean()) / expr.std(ddof=0)).values
    Model = RandomForestRegressor if method == "RF" else ExtraTreesRegressor

    W = np.zeros((p, p))
    for j in range(p):                       # j = target gene
        predictors = [i for i in range(p) if i != j]   # exclude self-regulation
        model = Model(n_estimators=n_trees, max_features=max_features,
                      random_state=random_state, n_jobs=n_jobs)
        model.fit(X_all[:, predictors], X_all[:, j])
        W[predictors, j] = _raw_importance(model, len(predictors))
    return pd.DataFrame(W, index=genes, columns=genes)


def pearson_baseline(expr):
    """Course-method baseline: |Pearson correlation| between genes.
    Symmetric, so it cannot tell i->j from j->i."""
    C = expr.corr().abs().to_numpy(copy=True)   # copy: pandas>=3 arrays are read-only
    np.fill_diagonal(C, 0.0)
    return pd.DataFrame(C, index=expr.columns, columns=expr.columns)
