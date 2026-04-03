import time
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, coo_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.stats import rankdata
from concurrent.futures import ProcessPoolExecutor
import networkx as nx
import CH_scores

ALL_MODELS = ['RA', 'CH1', 'CH2', 'CH3', 'iLCL', 'CH3.1']


def _replace_inf_distances(distance_matrix):
    distance_matrix = np.array(distance_matrix)
    inf_indices = np.isinf(distance_matrix)
    if not np.any(inf_indices):
        return distance_matrix
    G = nx.Graph(distance_matrix != np.inf)
    components = list(nx.connected_components(G))
    max_distances = []
    for component in components:
        idx = list(component)
        sub = distance_matrix[np.ix_(idx, idx)]
        max_distances.append(np.nanmax(sub))
    replaced = np.copy(distance_matrix)
    replaced[inf_indices] = np.sum(max_distances)
    return replaced


def _compute_SPcorr_and_rank_scores(x, S):
    distance = x * np.abs(S - np.min(S[S > 0]) - np.max(S))
    dist_matrix = shortest_path(csr_matrix(distance), method='D', directed=False)
    dist_matrix = _replace_inf_distances(dist_matrix)
    pcorr = np.corrcoef(dist_matrix)
    n = dist_matrix.shape[0]
    e2, e1 = np.tril_indices(n, -1)
    tmp = np.vstack((S[e1, e2], pcorr[e1, e2])).T
    _, unique_indices = np.unique(tmp, axis=0, return_inverse=True)
    unique_indices += 1
    tmp = rankdata(unique_indices, method='min')
    return coo_matrix((tmp, (e1, e2)), shape=(n, n))


def _rank_scores(x, S):
    # Rank by score only, without SPcorr — used for raw RA_L3
    distance = x * np.abs(S - np.min(S[S > 0]) - np.max(S))
    dist_matrix = shortest_path(csr_matrix(distance), method='D', directed=False)
    dist_matrix = _replace_inf_distances(dist_matrix)
    n = dist_matrix.shape[0]
    e2, e1 = np.tril_indices(n, -1)
    tmp = S[e1, e2]
    _, unique_indices = np.unique(tmp, return_inverse=True)
    unique_indices += 1
    tmp = rankdata(unique_indices, method='min')
    return coo_matrix((tmp, (e1, e2)), shape=(n, n))


def _subrank_task(args):
    method, x, S, rows, cols, is_ra_l3 = args
    print(f"starting subranking {method}...")
    t0 = time.time()
    if is_ra_l3:
        raw   = _rank_scores(x, S).toarray()[rows, cols]
        sprank = _compute_SPcorr_and_rank_scores(x, S).toarray()[rows, cols]
        print(f"  subranking {method}: {time.time() - t0:.2f}s")
        return [('RA_L3', raw), ('RA_L3_subranking', sprank)]
    else:
        result = _compute_SPcorr_and_rank_scores(x, S).toarray()[rows, cols]
        print(f"  subranking {method}: {time.time() - t0:.2f}s")
        return [(method, result)]


def cha_linkpred_monopartite(x, methods):
    if hasattr(x, 'toarray'):
        x = x.toarray()
    x = np.asarray(x, dtype=float)
    N = x.shape[0]

    lengths_set = set(int(m.rsplit('_L', 1)[1]) for m in methods)
    lengths = np.array(sorted(lengths_set), dtype=float)
    length_max = int(np.max(lengths))
    L = len(lengths)
    length_to_idx = {int(l): i for i, l in enumerate(lengths)}

    model_indexes = [float(i) for i in range(len(ALL_MODELS))]
    sparse_A = csr_matrix(x)
    scores_raw = CH_scores.CH_scores_new_v2(
        list(sparse_A.indices.astype(int)), list(sparse_A.indptr.astype(int)),
        N, list(lengths), L, length_max, model_indexes, len(model_indexes)
    )
    scores_raw = scores_raw.reshape((len(ALL_MODELS), L, N, N))

    cols, rows = np.tril_indices(N, -1)
    col_data = {'node1': rows + 1, 'node2': cols + 1}

    idx_ra_l3 = methods.index('RA_L3') if 'RA_L3' in methods else None

    tasks = []
    for i2 in range(L):
        path_len = int(lengths[i2])
        for i1, model_name in enumerate(ALL_MODELS):
            method = f'{model_name}_L{path_len}'
            if method not in methods:
                continue
            S = scores_raw[i1, i2, :, :]
            is_ra_l3 = (method == 'RA_L3' and idx_ra_l3 is not None)
            tasks.append((method, x, S, rows, cols, is_ra_l3))

    with ProcessPoolExecutor() as executor:
        for pairs in executor.map(_subrank_task, tasks):
            for col, data in pairs:
                col_data[col] = data

    return pd.DataFrame(col_data)
