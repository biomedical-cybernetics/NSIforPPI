% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script computes several topological measures 
% based on original network for Yeast DIP network.


networks = {'Yeast_DIP_net'};

measures = {'N', 'E', 'avgdeg', 'clustering', 'char_path', 'LCPcorr'};

for i = 1:length(networks)
    display(networks{i});
    load(['../matrix/' networks{i} '.mat'], 'x_lcc');
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
        st = compute_topological_measures(x_lcc, measures);
    end
    save(filename, 'st');
end
