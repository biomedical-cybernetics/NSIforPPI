% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script computes several topological measures 
% based on original network for various organisms.


networks = {'S_Pombe_BioGRID','Yeast','Mouse_BioGRID','Worm_BioGRID',...
    'HI_II_14','Arabidopsis','Yeast_BioGRID','Interactome3D','Lit_BM_13','Hein',...
    'Fly_BioGRID','Arabidopsis_BioGRID','HI_III','Bioplex','Human_BioGRID'};

% measures = {'N', 'E', 'avgdeg', 'density', 'clustering', ...
%     'char_path', 'efficiency_glob', 'efficiency_loc', 'closeness', ...
%     'EBC', 'BC', 'radiality', 'LCPcorr', 'assortativity', ...
%     'modularity', 'struct_cons', 'powerlaw_p', 'powerlaw_gamma', ...
%     'smallworld_omega', 'smallworld_sigma', ...
%     'smallworld_omega_eff', 'smallworld_sigma_eff'};

measures = {'N', 'E', 'avgdeg', 'density', 'clustering', ...
    'char_path', 'efficiency_glob', 'efficiency_loc', 'closeness', ...
    'EBC', 'BC', 'radiality', 'LCPcorr', 'assortativity', ...
    'struct_cons', 'powerlaw_p', 'powerlaw_gamma'};


for i = 1:length(networks)
    display(networks{i});
    if strcmp(networks{i}, 'Yeast_DIP_net')
        load(['../matrix/' networks{i} '.mat'], 'x_lcc');
    else
        load(['../matrix/' networks{i} '.mat'], 'x');
    end
    filename = ['../statistics/results/' networks{i} '_statistics.mat'];
    if exist(filename, 'file')
        load(filename, 'st');
        measures_diff = setdiff(measures, fieldnames(st));
        if isempty(measures_diff)
            continue
        end
        st2 = compute_topological_measures(x, measures_diff);
        for j = 1:length(measures_diff)
            eval(['st. ' measures_diff{j} ' = st2.' measures_diff{j} ';'])
        end
    else
        st = compute_topological_measures(x, measures);
    end
    save(filename, 'st');
end
