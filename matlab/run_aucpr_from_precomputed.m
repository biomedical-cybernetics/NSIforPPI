function run_aucpr_from_precomputed(net_file, perturbed_net_file, pairs_file, varargin)
% Run CH link prediction + AUCPR evaluation from pre-computed .mat files.
%
% Usage (default paths):
%   run_aucpr_from_precomputed()
%
% Usage (custom paths):
%   run_aucpr_from_precomputed('../matrix/Yeast_DIP_net.mat', ...
%       '../matrix/network_perturbed_10percent_GSP_GSN_Yeast_noconn.mat', ...
%       '../matrix/list_pairs_10percent_GSP_GSN_Yeast_noconn.mat')

if nargin < 1 || isempty(net_file)
    net_file = '../matrix/Yeast_DIP_net.mat';
end
if nargin < 2 || isempty(perturbed_net_file)
    perturbed_net_file = '../matrix/network_perturbed_10percent_GSP_GSN_Yeast_noconn.mat';
end
if nargin < 3 || isempty(pairs_file)
    pairs_file = '../matrix/list_pairs_10percent_GSP_GSN_Yeast_noconn.mat';
end

ip = inputParser;
addParameter(ip, 'output_file', '../table/table_auc_pr_CH_Yeast_DIP_net.xlsx');
parse(ip, varargin{:});
filename = ip.Results.output_file;

% ---- Load data ----
fprintf('Loading %s\n', net_file);
load(net_file, 'x_lcc');

fprintf('Loading %s\n', perturbed_net_file);
load(perturbed_net_file, 'L_net');

fprintf('Loading %s\n', pairs_file);
load(pairs_file, 'L_pairs');

num_sims = length(L_net);
fprintf('Loaded %d simulations.\n\n', num_sims);

% ---- CH link prediction ----
methods = {'RA_L2','CH1_L2','CH2_L2','CH3_L2','CH3.1_L2','iLCL_L2', ...
           'RA_L3','CH1_L3','CH2_L3','CH3_L3','CH3.1_L3','iLCL_L3'};

S = cell(num_sims, 1);
for i = 1:num_sims
    fprintf('CH link prediction: simulation %d/%d\n', i, num_sims);
    S{i} = CHA_linkpred_monopartite_final(L_net{i}, methods);
end

results_dir = '../network_similarities/CH_L2_L3/results/';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
save([results_dir 'CH_L2_L3_scores_net_perturbed.mat'], 'S', '-v7.3');
fprintf('CH scores saved.\n\n');

% ---- AUCPR evaluation ----
CH_methods = {'RA_L2','CH1_L2','CH2_L2','CH3_L2','CH3_1_L2','iLCL_L2', ...
              'RA_L3','RA_L3_subranking','CH1_L3','CH2_L3','CH3_L3','CH3_1_L3','iLCL_L3'};

res_CH_AUCPR = zeros(num_sims, length(CH_methods));

for i = 1:num_sims
    S_tmp = table2array(S{i});
    CH_mat = zeros(size(x_lcc));
    ind = sub2ind(size(CH_mat), S_tmp(:,1), S_tmp(:,2));
    P = L_pairs{i,1};
    N = L_pairs{i,2};
    labels = repelem([1, 0], size(P,1))';

    for j = 1:length(CH_methods)
        CH_mat(ind) = S_tmp(:, j+2);
        CH_mat_sym = CH_mat + CH_mat.';
        pos_scores = arrayfun(@(k) CH_mat_sym(P(k,1), P(k,2)), 1:size(P,1))';
        neg_scores = arrayfun(@(k) CH_mat_sym(N(k,1), N(k,2)), 1:size(N,1))';
        measures = prediction_evaluation([pos_scores; neg_scores], labels);
        res_CH_AUCPR(i, j) = measures.auc_pr;
    end
    fprintf('AUCPR evaluation: simulation %d/%d done\n', i, num_sims);
end

% ---- Build and write output table ----
tp = round(res_CH_AUCPR', 3);
tr = zeros(size(tp));
for j = 1:num_sims
    tr(:,j) = tiedrank(-tp(:,j));
end
t = [mean(tr,2), mean(tp,2), tp];
[t, idx] = sortrows(t, [1, -2]);

mean_ranking    = t(:,1)';
mean_aucpr      = t(:,2)';
sim_data        = t(:,3:end)';
col_methods     = CH_methods(idx);
sim_names       = arrayfun(@(k) sprintf('Realization_%02d', k), (1:num_sims)', 'UniformOutput', false);

table_dir = fileparts(filename);
if ~isempty(table_dir) && ~exist(table_dir, 'dir'), mkdir(table_dir); end

writetable(cell2table([{'methods'},      col_methods]),            filename, 'Range', 'A1', 'WriteVariableNames', false);
writetable(cell2table([{'mean_ranking'}, num2cell(mean_ranking)]), filename, 'Range', 'A2', 'WriteVariableNames', false);
writetable(cell2table([{'mean_aucpr'},   num2cell(mean_aucpr)]),   filename, 'Range', 'A3', 'WriteVariableNames', false);
writetable(cell2table([sim_names,        num2cell(sim_data)]),     filename, 'Range', 'A4', 'WriteVariableNames', false);

fprintf('\nResults written to %s\n', filename);
disp(round(mean(res_CH_AUCPR), 3));
end
