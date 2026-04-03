% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script generates the AUP evaluation table 
% of top-100 (AUP@100) interactions for PPI networks.

filename = '../table/table_aup_CHv2_PPI_top100.xlsx';

networks = {'S_Pombe_BioGRID','Yeast','Mouse_BioGRID','Worm_BioGRID',...
    'HI_II_14','Arabidopsis','Yeast_BioGRID','Interactome3D','Lit_BM_13','Hein',...
    'Fly_BioGRID','Arabidopsis_BioGRID','HI_III','Bioplex','Human_BioGRID'};

methods = {'RA_L2', 'CH1_L2', 'CH2_L2', 'CH3_L2', 'iLCL_L2', ...
           'RA_L3', 'CH1_L3', 'CH2_L3', 'CH3_L3', 'iLCL_L3'};

eval_metric = 'aup';
name_eval_metric = 'aup';
top = 100;


tp = zeros(length(methods),length(networks));
tr = zeros(length(methods),length(networks));
for j = 1:(length(networks))
    if strcmp(eval_metric, 'aup')
        for i = 1:length(methods)
            load(['../network_similarities/CH_L2_L3/results/' networks{j} '_linkrem10noconn_CHv2_L2_L3_linkpred.mat'], [methods{i} '_curve'])
            eval(['curve = ' methods{i} '_curve; clear ' methods{i} '_curve'])
            aup = zeros(length(curve),1);
            for k = 1:length(curve)
                aup(k) = trapz(1:top,curve{k}(1:top))/(top-1);
            end
            tp(i,j) = round(mean(aup),2);
        end
        tr(:,j) = tiedrank(-tp(:,j));
        
    elseif strcmp(eval_metric, 'auc_mroc')
        for i = 1:length(methods)
            load(['../network_similarities/CH_L2_L3/results/' networks{j} '_linkrem10noconn_CHv2_L2_L3_linkpred_AUC_mROC.mat'], [methods{i} '_auc_mroc'])
            eval(['tp(i,j) = ' 'round(mean('  methods{i} '_auc_mroc), 2);' ]);
        end
        tr(:,j) = tiedrank(-tp(:,j));
        
    elseif strcmp(eval_metric, 'auc_pr')
        for i = 1:length(methods)
            load(['../network_similarities/CH_L2_L3/results/' networks{j} '_linkrem10noconn_CHv2_L2_L3_linkpred_AUC_PR.mat'], [methods{i} '_auc_pr'])
            eval(['tp(i,j) = ' 'round(mean('  methods{i} '_auc_pr), 3);' ]);
        end     
        tr(:,j) = tiedrank(-tp(:,j));    
        
    end

end
    
t = [mean(tr,2) mean(tp,2) tp];
[t,idx] = sortrows(t,[1 -2]);
t = t';
rownames = [{'mean ranking'; ['mean ' name_eval_metric ]}; networks'];
writetable(cell2table(rownames), filename);
writetable(cell2table(methods(idx)), filename, 'Range', 'B1', 'WriteVariableNames', false);
writetable(array2table(t), filename, 'Range', 'B2', 'WriteVariableNames', false);