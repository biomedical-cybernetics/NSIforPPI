networks = {'S_Pombe_BioGRID','Yeast','Mouse_BioGRID','Worm_BioGRID',...
    'HI_II_14','Arabidopsis','Yeast_BioGRID','Interactome3D','Lit_BM_13','Hein',...
    'Fly_BioGRID','Arabidopsis_BioGRID','HI_III','Bioplex','Human_BioGRID'};

top = 100;
methods = {'RA_L3', 'CH1_L3' 'CH2_L3', 'CH3_L3'};
methods_plot = {'RA-L3', 'CH1-L3', 'CH2-L3', 'CH3-L3'};
nbars = length(methods);

pvalue = zeros(length(networks),1);
perc_improv = zeros(length(networks),1);
best = zeros(length(networks),1);
aupdiff = zeros(length(networks),1);
% This creates a cell array C with 15 rows and 1 column,
% where each cell contains an array of size 1x3 filled with zeros.
aup_mean_plot = repmat({zeros(1,3)}, 15, 1);

for j = 1:length(networks)    
    temp = load(['../pvalues_aup/' networks{j} '_pvalues_CHv1v2_aup_top100_fig2.mat'], 'p', 'methods', 'aup_mean');
    m1 = find(strcmp(methods{1}, temp.methods));
    m2 = find(strcmp(methods{2}, temp.methods));
    m3 = find(strcmp(methods{3}, temp.methods));
    m4 = find(strcmp(methods{4}, temp.methods));
    aup_mean_plot{j} = [temp.aup_mean(m1), temp.aup_mean(m2), temp.aup_mean(m3), temp.aup_mean(m4)];
    
    [~,imax] = max(temp.aup_mean([m1, m2, m3, m4]));
    [~,imin] = min(temp.aup_mean([m1, m2, m3, m4]));
    
    eval(['pvalue(j) = temp.p(m', num2str(imax), ',m1',');']);
    eval(['perc_improv(j) = round((temp.aup_mean(m', num2str(imax),') - temp.aup_mean(m1)) / temp.aup_mean(m1) * 100,1);']);
    best(j) = imax;
    eval(['aupdiff(j) = abs(temp.aup_mean(m', num2str(imax),') - temp.aup_mean(m', num2str(imin),'));']);

end

[~, idx] = sortrows([perc_improv pvalue], [-1 2]);
perc_improv = perc_improv(idx);
pvalue = pvalue(idx);
best = best(idx);
networks = networks(idx);
aup_mean_plot = aup_mean_plot(idx);

figure('color','white');
len = 5;
red = [1, 0, 0];
pink = [255, 192, 203]/255;
colors_p = [linspace(red(1),pink(1),len)', linspace(red(2),pink(2),len)', ...
 linspace(red(3),pink(3),len)'];
colors = {[0 0 0], colors_p(1,:), colors_p(2,:), colors_p(3,:), colors_p(4,:), colors_p(5,:)};
for j = 1:length(networks)
    subplot(5,3,j)
    aup = cell(length(methods),1);
    [y, idx] = sort(aup_mean_plot{j}, 'descend');
    b = bar(y, 'FaceColor','flat');
    m = methods(idx);
    xticklabels(methods_plot(idx));
    ylim([0 1])
    set(gca,'fontsize',13)
    set(gca, 'ytick', 0:0.2:1)
    yline(max(aup_mean_plot{j}), '--')
    ctr(1,:) = bsxfun(@plus, b(1).XData, b(1).XOffset');    
    ydt(1,:) = b(1).YData;                               
    text(ctr(1,:), ydt(1,:), sprintfc('%.3f', ydt(1,:)), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',14, 'Color','k')

    for i = 1:length(m)
        
        method = m{i};
        load(['../network_similarities/CH_L2_L3/results/' networks{j} '_linkrem10noconn_CHv2_L2_L3_linkpred.mat'], [method '_curve'])
        eval(['curve = ' method '_curve; clear ' method '_curve'])
        
        aup{i} = zeros(length(curve),1);
        for k = 1:length(curve)
            aup{i}(k) = trapz(1:top,curve{k}(1:top))/(top-1);
        end
        
        if ~strcmp(m{i}, 'RA_L3') && aup_mean_plot{j}(idx(i)) >= 0.8
            b.CData(i,:) = colors{2};
            if i == 1
                color_for_arrow = colors{2};
            end
        elseif ~strcmp(m{i}, 'RA_L3') && aup_mean_plot{j}(idx(i)) >= 0.6
            b.CData(i,:) = colors{3};
            if i == 1
                color_for_arrow = colors{3};
            end
        elseif ~strcmp(m{i}, 'RA_L3') && aup_mean_plot{j}(idx(i)) >= 0.4
            b.CData(i,:) = colors{4};
            if i == 1
                color_for_arrow = colors{4};
            end
        elseif ~strcmp(m{i}, 'RA_L3') && aup_mean_plot{j}(idx(i)) >= 0.2
            b.CData(i,:) = colors{5};
            if i == 1
                color_for_arrow = colors{5};
            end
        elseif ~strcmp(m{i}, 'RA_L3') && aup_mean_plot{j}(idx(i)) >= 0
            b.CData(i,:) = colors{6};
            if i == 1
                color_for_arrow = colors{6};
            end
        elseif strcmp(m{i}, 'RA_L3')
            b.CData(i,:) = colors{1};
            if i == 1
                color_for_arrow = colors{1};
            end
        end
        
        stderror = std(aup{i})/sqrt(size(aup{i}, 1));
        hold on
        err = errorbar(i, aup_mean_plot{j}(idx(i)), stderror, 'k', 'linestyle', 'none');
        err.LineWidth = 2.5;
        box off
        hold off
    
    end
    box on
    text(0.5, 1.1, strrep(networks{j},'_',' '), 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 15)
    text(0.80, 0.90, '\uparrow', 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 14, 'Color', color_for_arrow)
    if pvalue(j) <= 0.05
        text(0.78, 0.90, '*', 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 14, 'Color', color_for_arrow)
    end
    Pos = get(subplot(5,3,j),'Position');
    P1=[nbars, min(aup_mean_plot{j})];     
    P2=[nbars, max(aup_mean_plot{j})];         
    Xlim=[0 4.5];
    Ylim=[0 1];
    hold on
    X_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P1(1)-Xlim(1));
    X_conv(2)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P2(1)-Xlim(1));
    Y_conv(1)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P1(2)-Ylim(1));
    Y_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P2(2)-Ylim(1));
    xlim(Xlim)
    ylim(Ylim)
    annotation('arrow', X_conv, Y_conv, 'Color', [169 169 169]/255, 'LineWidth', 2)
   
    text(0.90, 0.90, sprintfc('%.1f%%', perc_improv(j)), 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 14, 'Color', color_for_arrow)
    if j==7
        ylabel('AUP', 'FontSize', 14)
    elseif j==14     
        xlabel('Methods', 'FontSize', 14)
    end
end
