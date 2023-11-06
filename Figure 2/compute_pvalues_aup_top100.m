% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script computes the permutation test for all the
% methods with subranking. 


networks = {'S_Pombe_BioGRID','Yeast','Mouse_BioGRID','Worm_BioGRID',...
    'HI_II_14','Arabidopsis','Yeast_BioGRID','Interactome3D','Lit_BM_13','Hein',...
    'Fly_BioGRID','Arabidopsis_BioGRID','HI_III','Bioplex','Human_BioGRID'};

methods = {'RA_L2', 'CH1_L2', 'CH2_L2', 'CH3_L2', ...
           'RA_L3', 'CH1_L3', 'CH2_L3', 'CH3_L3'};

top = 100;

for j = 1:length(networks)
    p = NaN(length(methods));
    aup_mean = zeros(length(methods),1);
    aup = cell(length(methods),1);
    for m = 1:length(methods)
        load(['../network_similarities/CH_L2_L3/results/' networks{j} '_linkrem10noconn_CHv2_L2_L3_linkpred.mat'], [methods{m} '_curve']);
        method = methods{m};
        eval(['curve = ' method '_curve; clear ' method '_curve'])
        aup{m} = zeros(length(curve),1);
        for k = 1:length(curve)
            aup{m}(k) = trapz(1:top,curve{k}(1:top))/(top-1);
        end
        aup_mean(m) = round(mean(aup{m}),2);
    end
    for m1 = 1:length(methods)-1
        for m2 = m1+1:length(methods)
            p(m1,m2) = permutation_test(aup{m1}, aup{m2}, 'mean', 'both', 10000, 1);
        end
    end
    p = max(p,p');
    save(['../pvalues_aup/' networks{j} '_pvalues_CHv1v2_aup_top100_fig2.mat'], 'p', 'aup_mean', 'methods')
end