% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates table of topological measures
% values based on original network for various organisms.


filename = '../statistics/statistics_networks_original';

networks = {'S_Pombe_BioGRID','Yeast','Mouse_BioGRID','Worm_BioGRID',...
    'HI_II_14','Arabidopsis','Yeast_BioGRID','Interactome3D','Lit_BM_13','Hein',...
    'Fly_BioGRID','Arabidopsis_BioGRID','HI_III','Bioplex','Human_BioGRID'};

measures = {'N', 'E', 'avgdeg', 'density', 'clustering', ...
    'char_path', 'efficiency_glob', 'efficiency_loc', 'closeness', ...
    'EBC', 'BC', 'radiality', 'LCPcorr', 'assortativity', ...
    'modularity', 'struct_cons', 'powerlaw_p', 'powerlaw_gamma'};

t = zeros(length(networks),length(measures));
for i = 1:length(networks)
    load(['../statistics/results/' networks{i} '_statistics.mat']);
    for j = 1:length(measures)
        eval(['t(i,j) = st.' measures{j} ';'])
    end
end

[~,idx] = sortrows(t, [1,2]);

xlswrite(filename, measures, 1, 'B1')
xlswrite(filename, networks(idx)', 1, 'A2')
xlswrite(filename, t(idx,:), 1, 'B2')
