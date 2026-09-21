#!/usr/bin/env python3
"""PAM50 subtype classifier from a 46-gene mutation panel.

Task: predict intrinsic subtype from non-silent mutation indicators.
Baselines -> logistic regression -> random forest, all evaluated with
stratified 5-fold cross-validation (cross_val_predict so every
prediction is out-of-fold). Writes metrics + figures to ml/results/.
"""

import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, confusion_matrix, f1_score,
                             precision_recall_fscore_support)
from sklearn.model_selection import StratifiedKFold, cross_val_predict

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "analysis" / "flagship-reproduction" / "data"
OUT = ROOT / "ml" / "results"
OUT.mkdir(parents=True, exist_ok=True)


def load_matrix():
    panel = json.loads((DATA / "panel_mutations.json").read_text())
    genes = panel["panel"]  # entrez -> symbol
    muts = panel["mutations"]  # entrez -> {sample_id: count}
    subs = pd.read_csv(DATA / "subtypes.tsv", sep="\t")
    patients = subs["patientId"].tolist()

    X = np.zeros((len(patients), len(genes)), dtype=int)
    pid_to_row = {p: i for i, p in enumerate(patients)}
    for col, (entrez, _sym) in enumerate(sorted(genes.items(),
                                                key=lambda kv: kv[1])):
        for sample_id in muts.get(entrez, {}):
            patient = sample_id[:12]  # TCGA barcode -> patient
            if patient in pid_to_row:
                X[pid_to_row[patient], col] = 1
    feat_names = [sym for _, sym in sorted(genes.items(), key=lambda kv: kv[1])]
    y = subs["pam50"].to_numpy()
    return X, y, feat_names, patients


def evaluate(name, model, X, y, cv):
    pred = cross_val_predict(model, X, y, cv=cv)
    acc = accuracy_score(y, pred)
    mf1 = f1_score(y, pred, average="macro")
    p, r, f, s = precision_recall_fscore_support(y, pred, zero_division=0)
    labels = sorted(set(y))
    return {
        "model": name, "accuracy": round(acc, 4), "macro_f1": round(mf1, 4),
        "per_class": {l: {"precision": round(float(pi), 3),
                          "recall": round(float(ri), 3),
                          "f1": round(float(fi), 3),
                          "support": int(si)}
                      for l, pi, ri, fi, si in zip(labels, p, r, f, s)},
        "confusion": confusion_matrix(y, pred, labels=labels).tolist(),
        "labels": labels,
        "_pred": pred,
    }


def main():
    X, y, feats, patients = load_matrix()
    print(f"matrix {X.shape}; class counts: {dict(Counter(y))}")

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    models = [
        ("dummy-most_frequent", DummyClassifier(strategy="most_frequent")),
        ("logreg-l2-balanced", LogisticRegression(
            C=1.0, class_weight="balanced", max_iter=2000)),
        ("rf-500-balanced", RandomForestClassifier(
            n_estimators=500, class_weight="balanced_subsample",
            random_state=42, n_jobs=-1)),
    ]
    results = [evaluate(n, m, X, y, cv) for n, m in models]
    for r in results:
        print(f"{r['model']}: acc={r['accuracy']} macroF1={r['macro_f1']}")

    # Full-fit RF for importances (reported as descriptive, not CV)
    rf = models[2][1].fit(X, y)
    imp = sorted(zip(feats, rf.feature_importances_),
                 key=lambda t: -t[1])[:15]

    # Learning curve: logistic macro-F1 vs train fraction (manual, CV-safe)
    lc = []
    for frac in (0.2, 0.4, 0.6, 0.8, 1.0):
        n = max(1, int(len(y) * frac))
        rng = np.random.RandomState(42)
        idx = rng.permutation(len(y))[:n]
        m = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)
        pred = cross_val_predict(m, X[idx], y[idx],
                                 cv=StratifiedKFold(5, shuffle=True,
                                                    random_state=42))
        lc.append({"train_n": n,
                   "macro_f1": round(f1_score(y[idx], pred, average="macro"), 4)})

    (OUT / "metrics.json").write_text(json.dumps({
        "n_patients": len(y), "n_features": X.shape[1],
        "feature_names": feats,
        "class_counts": {k: int(v) for k, v in Counter(y).items()},
        "models": [{k: v for k, v in r.items() if k != "_pred"}
                   for r in results],
        "rf_top15_importances": [(g, round(float(i), 4)) for g, i in imp],
        "learning_curve": lc,
        "cv": "StratifiedKFold(5, shuffle, seed=42), cross_val_predict",
    }, indent=1))

    # Figures
    best = max(results, key=lambda r: r["macro_f1"])
    labels = best["labels"]
    cm = np.array(best["confusion"])
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.imshow(cm, cmap="Blues")
    ax.set_xticks(range(len(labels)), [l.replace("BRCA_", "") for l in labels],
                  rotation=45, ha="right")
    ax.set_yticks(range(len(labels)), [l.replace("BRCA_", "") for l in labels])
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, cm[i, j], ha="center",
                    color="white" if cm[i, j] > cm.max() / 2 else "black")
    ax.set_xlabel("predicted"); ax.set_ylabel("true")
    ax.set_title(f"{best['model']} — out-of-fold confusion")
    fig.tight_layout(); fig.savefig(OUT / "confusion_matrix.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    names = [g for g, _ in imp][::-1]
    vals = [v for _, v in imp][::-1]
    ax.barh(names, vals); ax.set_xlabel("RF importance (full fit)")
    fig.tight_layout(); fig.savefig(OUT / "feature_importance.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot([p["train_n"] for p in lc], [p["macro_f1"] for p in lc], "o-")
    ax.set_xlabel("train n"); ax.set_ylabel("macro F1 (CV)")
    ax.set_title("LogReg learning curve")
    fig.tight_layout(); fig.savefig(OUT / "learning_curve.png", dpi=150)
    plt.close(fig)

    print(f"wrote {OUT}/metrics.json + 3 figures")


if __name__ == "__main__":
    main()
