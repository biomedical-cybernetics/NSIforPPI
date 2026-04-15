import argparse
import os
import pickle

import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import rankdata

from ch_linkpred import cha_linkpred_monopartite
from prediction_evaluation import prediction_evaluation

CH_METHODS = [
    'RA_L2', 'CH1_L2', 'CH2_L2', 'CH3_L2', 'iLCL_L2', 'CH3.1_L2',
    'RA_L3', 'RA_L3_subranking',
    'CH1_L3', 'CH2_L3', 'CH3_L3', 'iLCL_L3', 'CH3.1_L3',
]


def _load_mat(filepath):
    return loadmat(filepath)


def _load_field(data, candidates, label):
    for name in candidates:
        if name in data:
            return data[name]
    available = [k for k in data if not k.startswith('_')]
    raise KeyError(f"Variable not found in {label}. Available: {available}. Expected: {candidates}")


def _to_dense(mat):
    return mat.toarray() if hasattr(mat, 'toarray') else np.asarray(mat)


def run_aucpr_from_precomputed(
    net_file='../matrix/Yeast_DIP_net.mat',
    perturbed_file='../matrix/network_perturbed_10percent_GSP_GSN_Yeast_noconn.mat',
    pairs_file='../matrix/list_pairs_10percent_GSP_GSN_Yeast_noconn.mat',
    output_file='../table/table_auc_pr_CH_Yeast_DIP_net.xlsx',
    n_jobs=None,
):
    print(f"Loading network from {net_file}")
    net_data = _load_mat(net_file)
    x_lcc = _load_field(net_data, ['x_lcc', 'A', 'adj', 'network'], net_file)
    N = _to_dense(x_lcc).shape[0]
    print(f"  Network size: {N} nodes")

    print(f"Loading perturbed networks from {perturbed_file}")
    L_net = _load_field(_load_mat(perturbed_file), ['L_net', 'L_nets'], perturbed_file)

    print(f"Loading pairs from {pairs_file}")
    L_pairs = _load_field(_load_mat(pairs_file), ['L_pairs', 'pairs'], pairs_file)

    num_sims = L_net.shape[0]
    print(f"  {num_sims} simulations loaded\n")

    methods = ['RA_L2', 'CH1_L2', 'CH2_L2', 'CH3_L2', 'iLCL_L2', 'CH3.1_L2',
               'RA_L3', 'CH1_L3', 'CH2_L3', 'CH3_L3', 'iLCL_L3', 'CH3.1_L3']

    S_list = []
    for i in range(num_sims):
        print(f"CH link prediction: simulation {i+1}/{num_sims}")
        df = cha_linkpred_monopartite(_to_dense(L_net[i, 0]), methods, n_jobs=n_jobs)
        S_list.append(df)

    results_dir = '../network_similarities/CH_L2_L3/results/'
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, 'CH_L2_L3_scores_net_perturbed.pkl'), 'wb') as f:
        pickle.dump(S_list, f)
    print(f"CH scores saved.\n")

    res_CH_AUCPR = np.zeros((num_sims, len(CH_METHODS)))
    for i in range(num_sims):
        df = S_list[i]
        node1_0 = df['node1'].values - 1
        node2_0 = df['node2'].values - 1
        P_pairs = np.asarray(L_pairs[i, 0], dtype=int)
        N_pairs = np.asarray(L_pairs[i, 1], dtype=int)
        n_pos = P_pairs.shape[0]
        labels = np.concatenate([np.ones(n_pos), np.zeros(n_pos)])

        for j, method in enumerate(CH_METHODS):
            ch_mat = np.zeros((N, N))
            ch_mat[node1_0, node2_0] = df[method].values
            ch_mat_sym = ch_mat + ch_mat.T
            pos_scores = ch_mat_sym[P_pairs[:, 0] - 1, P_pairs[:, 1] - 1]
            neg_scores = ch_mat_sym[N_pairs[:, 0] - 1, N_pairs[:, 1] - 1]
            measures = prediction_evaluation(np.concatenate([pos_scores, neg_scores]), labels)
            res_CH_AUCPR[i, j] = measures['auc_pr']

        print(f"AUCPR evaluation: simulation {i+1}/{num_sims} done")

    tp = np.round(res_CH_AUCPR.T, 3)
    tr = np.zeros_like(tp)
    for j in range(num_sims):
        tr[:, j] = rankdata(-tp[:, j], method='average')
    t = np.column_stack([tr.mean(axis=1), tp.mean(axis=1), tp])
    sort_idx = np.lexsort((-t[:, 1], t[:, 0]))
    t = t[sort_idx]

    col_methods_sorted = [CH_METHODS[i] for i in sort_idx]
    sim_data = t[:, 2:].T
    sim_names = [f'Realization_{k+1:02d}' for k in range(num_sims)]

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    rows = (
        [['methods']      + col_methods_sorted] +
        [['mean_ranking'] + list(np.round(t[:, 0], 3))] +
        [['mean_aucpr']   + list(np.round(t[:, 1], 3))] +
        [[name] + list(np.round(row, 3)) for name, row in zip(sim_names, sim_data)]
    )
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        pd.DataFrame(rows).to_excel(writer, index=False, header=False)

    print(f"\nResults written to {output_file}")
    print(np.round(res_CH_AUCPR.mean(axis=0), 3))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--net_file',       default='../matrix/Yeast_DIP_net.mat')
    parser.add_argument('--perturbed_file', default='../matrix/network_perturbed_10percent_GSP_GSN_Yeast_noconn.mat')
    parser.add_argument('--pairs_file',     default='../matrix/list_pairs_10percent_GSP_GSN_Yeast_noconn.mat')
    parser.add_argument('--output_file',    default='../table/table_auc_pr_CH_Yeast_DIP_net.xlsx')
    parser.add_argument('--n_jobs',         type=int, default=None, help='Number of parallel workers (1=serial, None=all CPUs)')
    args = parser.parse_args()
    run_aucpr_from_precomputed(args.net_file, args.perturbed_file, args.pairs_file, args.output_file, args.n_jobs)
