% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script computes the CH similarity scores.


parpool();
load('../matrix/list_pairs_10percent_GSP_GSN_Yeast.mat');
load('../matrix/network_perturbed_10percent_GSP_GSN_Yeast.mat');
methods = {'RA_L2','CH1_L2','CH2_L2','CH3_L2','iLCL_L2','RA_L3','CH1_L3','CH2_L3','CH3_L3','iLCL_L3'};
CHA_option = {'CH2_L2', 'CH3_L2', 'CH2_L3', 'CH3_L3'};

S = cell(length(L_net), 1);
info = cell(length(L_net), 1);
for i = 1:length(L_net)
    [scores, CHA_info] = CHA_linkpred_monopartite(L_net{i}, methods, CHA_option);
    S{i} = scores;
    info{i} = CHA_info;
    
end

save('../network_similarities/CH_L2_L3/results/CHA_L3_scores_net_perturbed_10percent_GSP_GSN_Yeast_DIP.mat', 'S', 'info', '-v7.3');