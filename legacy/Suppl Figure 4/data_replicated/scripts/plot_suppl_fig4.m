% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script generates the precision curve evaluation
% for PPI networks: CH3-L3 vs RA-L3 methods.


networks = {'S_Pombe_BioGRID','Yeast','Mouse_BioGRID','Worm_BioGRID',...
    'HI_II_14','Arabidopsis','Yeast_BioGRID','Interactome3D','Lit_BM_13','Hein',...
    'Fly_BioGRID','Arabidopsis_BioGRID','HI_III','Bioplex','Human_BioGRID'};

top = 100;
methods = {'CH3_L3', 'RA_L3'};
methods_plot = {'CH3_L3', 'RA_L3'};
color = {'r', 'k'};
linetype = {'-', '-'};

pvalue = zeros(length(networks),1);
best = zeros(length(networks),1);
aupdiff = zeros(length(networks),1);
for j = 1:length(networks)
    temp = load(['../pvalues_aup/' networks{j} '_pvalues_CHv1v2_aup_top100_fig2.mat'], 'p', 'methods', 'aup_mean');
    m1 = find(strcmp(methods{1}, temp.methods));
    m2 = find(strcmp(methods{2}, temp.methods));
    pvalue(j) = temp.p(m1,m2);
    [~, best(j)] = max([temp.aup_mean(m1),temp.aup_mean(m2)]);
    aupdiff(j) = abs(temp.aup_mean(m1)-temp.aup_mean(m2));
 
end

[~, idx] = sortrows([aupdiff pvalue], [-1 2]);
pvalue = pvalue(idx);
best = best(idx);
networks = networks(idx);

figure('color','white')
for j = 1:length(networks)
    subplot(5,3,j)
    hold on
    for i = 1:length(methods)
        method = strrep(methods{i}, 'SP_v2', '');
        load(['../network_similarities/CH_L2_L3/results/' networks{j} '_linkrem10noconn_CHv2_L2_L3_linkpred.mat'], [method '_curve'])
        eval(['temp = ' method '_curve; clear ' method '_curve;'])

        curves = zeros(length(temp{1}),length(temp));
        for k = 1:length(temp)
            curves(:,k) = temp{k};
        end
        curve = mean(curves,2);
        plot(1:top,curve(1:top), linetype{i}, 'Color', color{i}, 'LineWidth', 2);
        ylim([0 1])
        set(gca, 'ytick', 0:0.2:1)
        set(gca, 'xtick', 0:25:100)
    end
    box on
    text(0.5, 1.1, strrep(networks{j},'_',' '), 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 12)
    text(0.95, 0.90, '\uparrow', 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 14, 'Color', color{best(j)})
    if pvalue(j) <= 0.05
        text(0.90, 0.90, '*', 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 14, 'Color', color{best(j)})
    end
    if j==7
        ylabel('Precision', 'FontSize', 12)
    elseif j==14
        xlabel('Links removed', 'FontSize', 12)
        ax_pos = get(gca,'Position');
        legend(strrep(methods_plot,'_','-'),'FontSize', 12, 'Position', [ax_pos(1) -0.02 ax_pos(3) ax_pos(4)], 'Orientation', 'horizontal')
    end
end
