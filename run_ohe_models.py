"""
run_ohe_models.py
-----------------
Runs the same family of classifiers used in prediction_tests.py, but on
in-memory arrays (one-hot encoded + centered sequences) instead of the
.scheme CSV files.

Put this file next to prediction_tests.py, then from your notebook:

    from run_ohe_models import run_models
    results = run_models(centered_sequences, label_for_training, prepro=3, threads=20)

It writes one hyper-parameter curve per algorithm (same style as the original
script), a combined model-comparison plot, and a results CSV.
"""

import os
import time

import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from sklearn import preprocessing, decomposition
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import (accuracy_score, f1_score, precision_score,
                             recall_score)

# reuse the reporting function already written in prediction_tests.py
from prediction_tests import metrics


# --------------------------------------------------------------------------
# data preparation
# --------------------------------------------------------------------------
def flatten_ohe(sequences):
    """(N, channels, length) one-hot -> (N, channels * length) float32."""
    seqs = np.asarray(sequences)
    if seqs.ndim != 3:
        raise ValueError(f"expected (N, channels, length), got {seqs.shape}")
    return seqs.reshape(seqs.shape[0], -1).astype(np.float32)


def preprocess(feature_vectors, prepro, n_components=100, seed=42):
    """Same options as prediction_tests.classification():
    1 = raw, 2 = scaling, 3 = PCA, 4 = scaling + PCA."""
    if prepro == 1:
        print("### Any")
        return feature_vectors

    if prepro == 2:
        print("### Scaling")
        scaler = preprocessing.StandardScaler().fit(feature_vectors)
        return scaler.transform(feature_vectors)

    if prepro in (3, 4):
        x = feature_vectors
        if prepro == 4:
            scaler = preprocessing.StandardScaler().fit(x)
            x = scaler.transform(x)
            scaler = None

        n_components = min(n_components, min(x.shape) - 1)
        pca = decomposition.PCA(n_components=n_components,
                                svd_solver='randomized',
                                random_state=seed)
        pca.fit(x)
        x = pca.transform(x)
        print("### PCA + Scaling" if prepro == 4 else "### PCA")
        print('X_PCA:', x.shape)
        print('explained variance kept: %.4f' % pca.explained_variance_ratio_.sum())
        pca = None
        return x

    raise ValueError("prepro must be 1, 2, 3 or 4")


# --------------------------------------------------------------------------
# generic hyper-parameter sweep (replaces the 8 copy-pasted blocks)
# --------------------------------------------------------------------------
def _sweep(tag, label, build, values, xlabel, data, outdir, logx=False):
    X_train, Y_train, X_validation, Y_validation = data

    ytrain, yValidation = [], []
    for v in values:
        t0 = time.time()
        model = build(v)
        model.fit(X_train, Y_train)
        ytrain.append(f1_score(Y_train, model.predict(X_train), average='macro'))
        yValidation.append(f1_score(Y_validation, model.predict(X_validation),
                                    average='macro'))
        print(f'{tag} {xlabel}={v}  train={ytrain[-1]:.4f}  '
              f'validation={yValidation[-1]:.4f}  ({time.time() - t0:.1f}s)')

    plt.close('all')
    plt.figure(figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')
    plot = plt.semilogx if logx else plt.plot
    plot(values, ytrain, '-', label='Train')
    plot(values, yValidation, '-', label='Validation')
    plt.ylim((0, 1.1))
    plt.xlabel(xlabel)
    plt.ylabel('F1-Score')
    plt.legend()
    plt.savefig(os.path.join(outdir, f'{tag}-Algorithm_ohe.png'), dpi=100)
    plt.close('all')

    best_i = int(np.argmax(yValidation))
    best_v = values[best_i]
    print(f'### {tag} ### The best score with data validation: '
          f'{yValidation[best_i]}  with {xlabel}: {best_v}')

    best = build(best_v)
    best.fit(X_train, Y_train)
    predictions = best.predict(X_validation)
    metrics(Y_validation, predictions)          # from prediction_tests.py

    return {
        'algorithm': tag,
        'name': label,
        'best_param': f'{xlabel}={best_v}',
        'f1': f1_score(Y_validation, predictions, average='macro'),
        'accuracy': accuracy_score(Y_validation, predictions),
        'recall': recall_score(Y_validation, predictions, average='macro'),
        'precision': precision_score(Y_validation, predictions, average='macro'),
    }


# --------------------------------------------------------------------------
# comparison plot
# --------------------------------------------------------------------------
def comparison_plot(results, outdir, title='Model comparison (validation set)'):
    df = pd.DataFrame(results).sort_values('f1', ascending=True)

    f1_color, acc_color = '#3B6FD4', '#E08A1E'
    y = np.arange(len(df))
    h = 0.38

    plt.close('all')
    path = os.path.join(outdir, 'MODEL-comparison_ohe.png')

    # rc_context so seaborn's global style (set in prediction_tests.py)
    # does not repaint this figure
    with mpl.rc_context(mpl.rcParamsDefault):
        fig, ax = plt.subplots(figsize=(11, 0.95 * len(df) + 2.2), dpi=100)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

        ax.barh(y + h / 2, df['f1'], height=h, color=f1_color,
                label='F1 (macro)', zorder=3)
        ax.barh(y - h / 2, df['accuracy'], height=h, color=acc_color,
                label='Accuracy', zorder=3)

        # direct labels, so the bars never have to be read off the axis
        for yi, (f1v, accv) in enumerate(zip(df['f1'], df['accuracy'])):
            ax.text(f1v + 0.012, yi + h / 2, f'{f1v:.3f}', va='center',
                    fontsize=9, color='#3d3d3d')
            ax.text(accv + 0.012, yi - h / 2, f'{accv:.3f}', va='center',
                    fontsize=9, color='#3d3d3d')

        ax.set_yticks(y)
        ax.set_yticklabels(df['name'])
        ax.set_xlim(0, 1.12)
        ax.set_xlabel('Score')
        ax.set_title(title, loc='left', fontsize=13, pad=14)
        ax.xaxis.grid(True, color='#e3e3e3', linewidth=0.8, zorder=0)
        ax.yaxis.grid(False)
        ax.set_axisbelow(True)
        for side in ('top', 'right', 'left'):
            ax.spines[side].set_visible(False)
        ax.spines['bottom'].set_color('#cfcfcf')
        ax.tick_params(length=0)
        ax.legend(frameon=False, loc='lower right', ncol=2)

        fig.tight_layout()
        fig.savefig(path, dpi=130, facecolor='white')
    plt.close('all')
    print(f'comparison plot saved to {path}')
    return path


# --------------------------------------------------------------------------
# main entry point
# --------------------------------------------------------------------------
def run_models(sequences, labels, prepro=3, threads=20, n_components=100,
               validation_size=0.2, seed=7, outdir='.',
               algorithms=('LR', 'LDA', 'KNN', 'MLP', 'RF', 'DT', 'NB', 'SVC')):
    """sequences: (N, channels, length) one-hot array (centered_sequences)
       labels:    list/array of 0/1 labels (label_for_training)"""
    os.makedirs(outdir, exist_ok=True)

    feature_vectors = flatten_ohe(sequences)
    label_vectors = np.asarray(labels)
    print('flattened features:', feature_vectors.shape)
    print('labels:', np.bincount(label_vectors))

    x_data = preprocess(feature_vectors, prepro, n_components=n_components)
    feature_vectors = None

    X_train, X_validation, Y_train, Y_validation = train_test_split(
        x_data, list(label_vectors), test_size=validation_size,
        random_state=seed, stratify=list(label_vectors))
    x_data = None
    print('train:', X_train.shape, ' validation:', X_validation.shape)

    data = (X_train, Y_train, X_validation, Y_validation)
    results = []

    if 'LR' in algorithms:
        results.append(_sweep(
            'LR', 'Logistic Regression',
            lambda c: LogisticRegression(C=c, n_jobs=threads, max_iter=1000),
            [round(c, 2) for c in np.arange(0.1, 1.0, 0.1)], 'C', data, outdir))

    if 'LDA' in algorithms:
        results.append(_sweep(
            'LDA', 'Linear Discriminant Analysis',
            lambda t: LinearDiscriminantAnalysis(tol=t),
            [round(t, 5) for t in np.arange(1e-4, 1.1e-3, 1e-4)], 'tol',
            data, outdir))

    if 'KNN' in algorithms:
        results.append(_sweep(
            'KNN', 'K-Nearest Neighbors',
            lambda k: KNeighborsClassifier(n_neighbors=k, n_jobs=threads),
            list(range(1, 100, 10)), 'n-Neighbors', data, outdir))

    if 'MLP' in algorithms:
        results.append(_sweep(
            'MLP', 'Multilayer Perceptron',
            lambda n: MLPClassifier(solver='lbfgs', alpha=.5,
                                    hidden_layer_sizes=(n,), random_state=seed),
            list(range(50, 500, 50)), 'Neurons', data, outdir))

    if 'RF' in algorithms:
        results.append(_sweep(
            'RF', 'Random Forest',
            lambda n: RandomForestClassifier(n_estimators=n, n_jobs=threads,
                                             random_state=seed),
            list(range(10, 100, 10)), 'Trees', data, outdir))

    if 'DT' in algorithms:
        results.append(_sweep(
            'DT', 'Decision Tree',
            lambda d: DecisionTreeClassifier(max_depth=d, random_state=seed),
            list(range(1, 10)), 'Max Depth', data, outdir))

    if 'NB' in algorithms:
        results.append(_sweep(
            'NB', 'Gaussian Naive Bayes',
            lambda v: GaussianNB(var_smoothing=v),
            [1e-1, 1e-3, 1e-5, 1e-7, 1e-9, 1e-11, 1e-13, 1e-15, 1e-17, 1e-19],
            'var_smoothing', data, outdir, logx=True))

    if 'SVC' in algorithms:
        results.append(_sweep(
            'SVC', 'Support Vector Machine',
            lambda c: OneVsRestClassifier(SVC(C=c, gamma=1e-6), n_jobs=threads),
            list(range(10, 100, 10)), 'C', data, outdir))

    df = pd.DataFrame(results).sort_values('f1', ascending=False)
    csv_path = os.path.join(outdir, 'model_comparison_ohe.csv')
    df.to_csv(csv_path, index=False)
    print('\n### SUMMARY ###')
    print(df.to_string(index=False))
    print(f'results saved to {csv_path}')

    comparison_plot(results, outdir)
    return df


if __name__ == '__main__':
    raise SystemExit("import run_models() from your notebook instead of "
                     "running this file directly")