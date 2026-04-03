import numpy as np
from scipy.stats import rankdata


def prediction_evaluation(scores, labels):
    scores = np.asarray(scores, dtype=float).ravel()
    labels = np.asarray(labels, dtype=float).ravel()
    S = len(scores)
    P = int(np.sum(labels == 1))
    N = S - P
    if P == 0 or N == 0:
        raise ValueError("labels cannot be all ones or all zeros")
    measures = _compute_curves_measures(scores, labels, S, P, N)
    measures['ndcg'] = _compute_ndcg(scores, labels, P)
    measures['mcc'] = _compute_mcc(scores, labels, P, N)
    return measures


def _compute_curves_measures(scores, labels, S, P, N):
    idx = np.argsort(-scores, kind='stable')
    neg_scores_sorted = -scores[idx]
    labels_sorted = labels[idx]

    _, ut_first = np.unique(neg_scores_sorted, return_index=True)
    ut = np.concatenate([ut_first[1:] - 1, [S - 1]])

    tp = np.cumsum(labels_sorted)
    fp = np.cumsum(1.0 - labels_sorted)
    tp_rand = fp * (P / N)

    rng = np.arange(1, S + 1, dtype=float)
    prec = tp / rng
    tpr = tp / P
    fpr = fp / N
    tpr_m = np.log1p(tp) / np.log1p(P)
    fpr_m = np.log1p(fp) / np.log1p(N)
    tpr_m_rand = np.log1p(tp_rand) / np.log1p(P)
    denom = 1.0 - tpr_m_rand
    tpr_m_norm = np.where(denom == 0, 1.0,
                          (tpr_m - tpr_m_rand) / denom * (1.0 - fpr_m) + fpr_m)

    measures = {}
    measures['prec'] = float(prec[P - 1])
    if P == 1:
        measures['auc_prec'] = float(prec[0])
    else:
        measures['auc_prec'] = float(np.trapz(prec[:P], np.arange(1, P + 1)) / (P - 1))

    prec_ut = prec[ut]
    tpr_ut = np.concatenate([[0.0], tpr[ut]])
    fpr_ut = np.concatenate([[0.0], fpr[ut]])
    fpr_m_ut = np.concatenate([[0.0], fpr_m[ut]])
    tpr_m_norm_ut = np.concatenate([[0.0], tpr_m_norm[ut]])

    if np.all(tpr_ut[1:] == 1.0):
        measures['auc_pr'] = float(prec_ut[0])
    else:
        measures['auc_pr'] = float(np.trapz(prec_ut, tpr_ut[1:]) / (1.0 - tpr_ut[1]))

    measures['auc_roc'] = float(np.trapz(tpr_ut, fpr_ut))
    measures['auc_mroc'] = float(np.trapz(tpr_m_norm_ut, fpr_m_ut))
    return measures


def _compute_ndcg(scores, labels, P):
    ranks = rankdata(-scores, method='average')
    pos_ranks = ranks[labels == 1]
    dcg = float(np.sum(1.0 / np.log2(1.0 + pos_ranks)))
    idcg = float(np.sum(1.0 / np.log2(1.0 + np.arange(1, P + 1))))
    return dcg / idcg


def _compute_mcc(scores, labels, P, N):
    idx = np.argsort(-scores, kind='stable')
    labels_sorted = labels[idx]
    tp = int(np.sum(labels_sorted[:P]))
    fp = P - tp
    tn = N - fp
    fn = fp
    denom = float(np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    if denom == 0:
        return 0.0
    mcc = (tp * tn - fp * fn) / denom
    return 0.0 if (np.isinf(mcc) or np.isnan(mcc)) else float(mcc)
