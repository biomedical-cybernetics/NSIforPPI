% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script calculates the relative deviation of
% topological measures for PPI networks.


original_stats = readtable('../statistics/statistics_networks_original.xls');
linkrem10_stats = readtable('../statistics/statistics_networks_linkrem10.xls');
linkrem50_stats = readtable('../statistics/statistics_networks_linkrem50.xls');

original_stats = table2cell(original_stats);
linkrem10_stats = table2cell(linkrem10_stats);
linkrem50_stats = table2cell(linkrem50_stats);

networks = {'S_Pombe_BioGRID','Yeast','Mouse_BioGRID','Worm_BioGRID',...
    'HI_II_14','Arabidopsis','Yeast_BioGRID','Interactome3D','Lit_BM_13','Hein',...
    'Fly_BioGRID','Arabidopsis_BioGRID','HI_III','Bioplex','Human_BioGRID'};

measures = {'avgdeg', 'density', 'clustering', ...
    'char_path', 'efficiency_glob', 'efficiency_loc', 'closeness', ...
    'EBC', 'BC', 'radiality', 'LCPcorr', 'assortativity', ...
    'modularity', 'struct_cons', 'powerlaw_p', 'powerlaw_gamma'};

measures_to_plot = {'avgdeg', 'clustering', 'LCPcorr', 'char_path', ...
    'efficiency_glob', 'closeness', 'assortativity', 'powerlaw_gamma', ...
    'struct_cons', 'EBC', 'BC', 'modularity'};

measure_name_in_plot = {'average degree', 'clustering', 'LCPcorr', 'characteristic path length', ...
    'efficiency', 'closeness centrality', 'assortativity', 'powerlaw exponant', ...
    'structural consistency', 'EBC', 'BC', 'modularity'};

% Use cellfun to find indices of elements in A in B
indices = cell2mat(cellfun(@(x) find(ismember(measures, x)), measures_to_plot, 'UniformOutput', false));
original_stats = original_stats(:, indices+3);
linkrem10_stats = linkrem10_stats(:, indices+3);
linkrem50_stats = linkrem50_stats(:, indices+3);

dev_from_ori = zeros(length(measures), 2);
stderr_dev_from_ori = zeros(length(measures), 2);
for i = 1:length(measures_to_plot)
    dev_from_ori(i, 1) = round(mean(abs( ([original_stats{:,i}] - [linkrem10_stats{:,i}]) ./ [original_stats{:,i}])) * 100, 1);
    dev_from_ori(i, 2) = round(mean(abs( ([original_stats{:,i}] - [linkrem50_stats{:,i}]) ./ [original_stats{:,i}])) * 100, 1);
    
    stderr_dev_from_ori(i, 1) = round(std(abs(([original_stats{:,i}] - [linkrem10_stats{:,i}]) ./ [original_stats{:,i}])) * 100 /sqrt(size(original_stats,1)), 3);
    stderr_dev_from_ori(i, 2) = round(std(abs(([original_stats{:,i}] - [linkrem50_stats{:,i}]) ./ [original_stats{:,i}])) * 100 /sqrt(size(original_stats,1)), 3);
    
    subplot(4,3,i)
    b = bar(dev_from_ori(i,:), 'FaceColor', 'flat');
    b.CData(1,:) = [0 0 0];
    b.CData(2,:) = [1 0 0];
    title(measure_name_in_plot{i}, 'FontWeight','Normal');
    set(gca, 'xticklabel', {'10% removal', '50% removal'}, 'FontSize', 15);
    set(gca,'yticklabel',[]);
    ylabel('Deviation from original');
    ctr(1,:) = bsxfun(@plus, b(1).XData, b(1).XOffset');
    ydt(1,:) = b(1).YData;
    text(ctr(1,:), ydt(1,:) - 0.06, {[num2str(dev_from_ori(i,1)), '%'] ;[num2str(dev_from_ori(i,2)), '%' ] }, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15);
    hold on
    err =  errorbar(1:2, dev_from_ori(i, :), stderr_dev_from_ori(i, :), 'k', 'linestyle', 'none');
    err.LineWidth = 2.5;
    hold off
    
 
end