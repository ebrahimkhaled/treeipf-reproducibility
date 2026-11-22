# train scikit-learn random forests on every csv in a folder and save the trees as json
# (the json keeps only what a shared model file holds: splits, leaf values, node counts)
#
# usage: python sklearn_export.py <csv folder> <json folder>
# each csv has columns X1..Xp and y; file names start with "gaussian" or "binomial"

import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor


def tree_dict(est, task):
    t = est.tree_
    if task == "classification":
        v = t.value[:, 0, :]
        v = v / v.sum(axis=1, keepdims=True)
        value = v[:, 1]
    else:
        value = t.value[:, 0, 0]
    return {
        "feature": t.feature.tolist(),
        "threshold": t.threshold.tolist(),
        "left": t.children_left.tolist(),
        "right": t.children_right.tolist(),
        "value": value.tolist(),
        "n": t.n_node_samples.tolist(),
        "wn": t.weighted_n_node_samples.tolist(),
        "impurity": t.impurity.tolist(),
    }


def export_forest(csv_path, out_path, seed):
    d = pd.read_csv(csv_path)
    X = d.drop(columns="y")
    y = d["y"].to_numpy()
    task = "classification" if os.path.basename(csv_path).startswith("binomial") else "regression"
    common = dict(n_estimators=200, max_features=0.5, min_samples_leaf=10,
                  bootstrap=True, random_state=seed, n_jobs=1)
    if task == "classification":
        rf = RandomForestClassifier(**common).fit(X, y)
    else:
        rf = RandomForestRegressor(**common).fit(X, y)
    out = {
        "task": task,
        "feature_names": list(X.columns),
        "max_features": 0.5,
        "max_samples": None,
        "min_samples_leaf": 10,
        "trees": [tree_dict(e, task) for e in rf.estimators_],
    }
    with open(out_path, "w") as f:
        json.dump(out, f)


if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    os.makedirs(dst, exist_ok=True)
    files = sorted(f for f in os.listdir(src) if f.endswith(".csv"))
    for k, f in enumerate(files):
        out = os.path.join(dst, f.replace(".csv", ".json"))
        if os.path.exists(out):
            continue
        export_forest(os.path.join(src, f), out, seed=1000 + k)
    print("exported", len(files), "forests")
