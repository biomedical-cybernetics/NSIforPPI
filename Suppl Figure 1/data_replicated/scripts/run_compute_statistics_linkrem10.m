% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script computes several topological measures 
% based on perturbed network (10% of links deleted) for various organisms.


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
    load(['../sparsified_matrices/' networks{i} '_linkrem10.mat'], 'matrices');

    filename = ['../statistics/results/' networks{i} '_linkrem10_statistics.mat'];
    if exist(filename, 'file')
        load(filename, 'st');
        measures_diff = setdiff(measures, fieldnames(st{1}));
        if isempty(measures_diff)
            continue
        end
        st2 = cell(length(matrices),1);
        parfor m = 1:length(matrices)
            st2{m} = compute_topological_measures(matrices{m}, measures_diff);
        end
        for m = 1:length(matrices)
            for j = 1:length(measures_diff)
                eval(['st{m}. ' measures_diff{j} ' = st2{m}.' measures_diff{j} ';'])
            end
        end
    else
        st = cell(length(matrices),1);
        parfor m = 1:length(matrices)
            st{m} = compute_topological_measures(matrices{m}, measures);
        end
    end
    save(filename, 'st');
end
