% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script computes the CH similarity scores.


load('../matrix/Yeast_DIP_net.mat');
methods = {'RA_L2','CH1_L2','CH2_L2','CH3_L2','iLCL_L2','RA_L3','CH1_L3','CH2_L3','CH3_L3','iLCL_L3'};
CHA_option = {'CH2_L2', 'CH3_L2', 'CH2_L3', 'CH3_L3'};

[scores, CHA_info] = CHA_linkpred_monopartite(x_lcc, methods, CHA_option);

save('../network_similarities/CH_L2_L3/results/Yeast_DIP_CHv2_L2_L3.mat', 'scores', 'CHA_info', '-v7.3');