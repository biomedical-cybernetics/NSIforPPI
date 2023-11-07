% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script generates the table of network statistics for the PPI
% networks


filename = '../statistics/table_statistics_PPI.xlsx';

networks = {'S_Pombe_BioGRID','Yeast','Mouse_BioGRID','Worm_BioGRID',...
    'HI_II_14','Arabidopsis','Yeast_BioGRID','Interactome3D','Lit_BM_13','Hein',...
    'Fly_BioGRID','Arabidopsis_BioGRID','HI_III','Bioplex','Human_BioGRID','Human_Menche_2015'};

measures = {'N', 'E', 'avgdeg', 'clustering', 'char_path', 'LCPcorr'};

t = zeros(length(networks),length(measures));
for i = 1:length(networks)
    load(['../statistics/results/' networks{i} '_statistics.mat']);
    for j = 1:length(measures)
         eval(['t(i,j) = st.' measures{j} ';']) 
    end
end

[~,idx] = sortrows(t, [1,2]);

writetable(cell2table(networks(idx)'), filename, 'Range', 'A2', 'WriteVariableNames', false);
writetable(cell2table(measures), filename, 'Range', 'B1', 'WriteVariableNames', false);
writetable(array2table(t(idx,:)), filename, 'Range', 'B2', 'WriteVariableNames', false);