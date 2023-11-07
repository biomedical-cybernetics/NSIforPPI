% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script generates the AUC-mROC evaluation of 
% all non-observed interactions for PPI networks.


networks = {'S_Pombe_BioGRID','Yeast','Mouse_BioGRID','Worm_BioGRID',...
    'HI_II_14','Arabidopsis','Yeast_BioGRID','Interactome3D','Lit_BM_13','Hein',...
    'Fly_BioGRID','Arabidopsis_BioGRID','HI_III','Bioplex','Human_BioGRID'};

methods = {'RA_L2','CH1_L2','CH2_L2','CH3_L2', 'iLCL_L2', 'RA_L3', 'CH1_L3',...
    'CH2_L3','CH3_L3', 'iLCL_L3'};
CHA_option = {'CH2_L2', 'CH3_L2', 'CH2_L3', 'CH3_L3'};
iters = 10;

for i = 1:length(networks)
    
    % check output
    temp = [networks{i} '_linkrem10noconn_CHv2_L2_L3_linkpred_AUCmROC.TEMP'];
    output = ['../network_similarities/CH_L2_L3/results/' networks{i} '_linkrem10noconn_CHv2_L2_L3_linkpred_AUC_mROC.mat'];
    if exist(temp,'file') || exist(output,'file')
        continue;
    end
    fclose(fopen(temp, 'w'));
    
    % initialization
    load(['../matrix/' networks{i} '.mat'], 'x'); x = x;
    load(['../sparsified_matrices/' networks{i} '_linkrem10.mat'], 'matrices')
    fprintf('%d/%d: %s (N=%d) ... ', i, length(networks), networks{i}, length(x))
    vartosave = {};
    for m = 1:length(methods)
        eval([methods{m} '_auc_mroc = zeros(iters,1);'])
        vartosave = [vartosave {[methods{m} '_auc_mroc']}];
    end
    check_output = ['../../../Figure 2/data_replicated/network_similarities/CH_L2_L3/results/', networks{i}, '_similarity_scores.mat'];
    if ~exist(check_output,'file')
        % link prediction
        for j = 1:iters
            fprintf('%d ', j)
            
            [scores, CHA_info] = CHA_linkpred_monopartite(matrices{j}, methods, CHA_option);
            for m = 1:length(methods)
                method = methods{m};
                ranking = table2array(sortrows(scores(:, [1:2 3+m]), -3));
                % Indices of missing links
                [e1,e2] = find(triu(matrices{j}==0,1));
                % Find common rows
                [common_rows, ~, idx] = intersect([e1 e2], [ranking(:, 1) ranking(:, 2)], 'rows');
                % Get indices in ranking arrays
                idx_A = find(ismember([ranking(:, 1) ranking(:, 2)], common_rows, 'rows'));
                eval([method '_auc_mroc(j) = linkpred_evaluation(x, matrices{j}, ranking(idx_A,:));'])
            end
            clear scores CHA_info;
        end
        
        save(output, vartosave{:}, '-v7.3')
        delete(temp);
        
    else
        copyfile(check_output, '../network_similarities/CH_L2_L3/results/');
        load(['../network_similarities/CH_L2_L3/results/', networks{i}, '_similarity_scores.mat'], 'S');
        for j = 1:iters
            fprintf('%d ', j)
            scores = S{j};
            
            for m = 1:length(methods)
                method = methods{m};
                ranking = table2array(sortrows(scores(:, [1:2 3+m]), -3));
                % Indices of missing links
                [e1,e2] = find(triu(matrices{j}==0,1));
                % Find common rows
                [common_rows, ~, idx] = intersect([e1 e2], [ranking(:, 1) ranking(:, 2)], 'rows');
                % Get indices in ranking arrays
                idx_A = find(ismember([ranking(:, 1) ranking(:, 2)], common_rows, 'rows'));
                eval([method '_auc_mroc(j) = linkpred_evaluation(x, matrices{j}, ranking(idx_A,:));'])
            end
            clear scores;
            
        end
        save(output, vartosave{:}, '-v7.3')
        delete(temp);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [auc_mroc] = linkpred_evaluation(x, x_spars, ranking)

%%% INPUT %%%
% x       - original matrix
% x_spars - sparsified matrix
% ranking - ranking of the non-observed links

%%% OUTPUT %%%
% auc_mroc - area under curve mROC


% links removed
removed = x - x_spars;

% for each ranked link, 1 if correct and 0 if not
result = removed(sub2ind(size(removed),ranking(:,1),ranking(:,2)));

% compute the precision curve
[measures_CH,~,~] = prediction_evaluation(ranking(:,3), result);
auc_mroc = measures_CH.auc_mroc;

end
