% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script For each simulation (sparsified network)
% of a network, computes the CH similarity scores for all the methods.
% Also, computes the AUP evaluation for all the methods. 


networks = {'S_Pombe_BioGRID','Yeast','Mouse_BioGRID','Worm_BioGRID',...
    'HI_II_14','Arabidopsis','Yeast_BioGRID','Interactome3D','Lit_BM_13','Hein',...
    'Fly_BioGRID','Arabidopsis_BioGRID','HI_III','Bioplex','Human_BioGRID'};

methods = {'RA_L2', 'CH1_L2', 'CH2_L2', 'CH3_L2', 'iLCL_L2', 'RA_L3', 'CH1_L3', 'CH2_L3', 'CH3_L3', 'iLCL_L3'};
CHA_option = {'CH2_L2', 'CH3_L2', 'CH2_L3', 'CH3_L3'};
iters = 10;

for i = 1:length(networks)
    
    S = cell(iters,1);
    
    % check output
    temp = [networks{i} '_linkrem10noconn_CHv2_L2_L3_linkpred.TEMP'];
    output = ['../network_similarities/CH_L2_L3/results/' networks{i} '_linkrem10noconn_CHv2_L2_L3_linkpred.mat'];
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
        eval([methods{m} '_curve = cell(iters,1);'])
        eval([methods{m} '_p = zeros(iters,1);'])
        eval([methods{m} '_aup = zeros(iters,1);'])
        vartosave = [vartosave {[methods{m} '_curve'],[methods{m} '_p'],[methods{m} '_aup']}];
    end

    % link prediction
    for j = 1:iters
		fprintf('%d ', j)
        [scores, CHA_info] = CHA_linkpred_monopartite(matrices{j}, methods, CHA_option);
        S{j} = scores;
        for m = 1:length(methods)
            method = methods{m};
            ranking = table2array(sortrows(scores(:, [1:2 3+m]), -3));
            % Indices of missing links
            [e1,e2] = find(triu(matrices{j}==0,1));
            % Find common rows
            [common_rows, ~, idx] = intersect([e1 e2], [ranking(:, 1) ranking(:, 2)], 'rows');
            % Get indices in ranking arrays
            idx_A = find(ismember([ranking(:, 1) ranking(:, 2)], common_rows, 'rows'));
            eval(['[~, ' method '_curve{j}, ' method '_p(j), ' method '_aup(j)] = linkpred_evaluation(x, matrices{j}, ranking(idx_A,:));'])
        end
        clear scores CHA_info;
    end        
    save(['../network_similarities/CH_L2_L3/results/' networks{i}, '_similarity_scores.mat'], 'S', '-v7.3')
    save(output, vartosave{:}, '-v7.3')
    delete(temp);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ranking, curve, p, aup] = linkpred_evaluation(x, x_spars, ranking)

%%% INPUT %%%
% x       - original matrix
% x_spars - sparsified matrix
% ranking - ranking of the non-observed links

%%% OUTPUT %%%
% ranking - top-ranking
% curve - precision curve
% p - precision
% aup - area under precision

% links removed
removed = x - x_spars;
m = nnz(triu(removed,1));

% top-ranking
ranking = ranking(1:m,:);

% for each top-ranked link, 1 if correct and 0 if not
result = removed(sub2ind(size(removed),ranking(:,1),ranking(:,2)));

% compute the precision curve
curve = full(cumsum(result) ./ (1:m)');
p = curve(end);
aup = trapz(1:m,curve) / (m-1);

end
