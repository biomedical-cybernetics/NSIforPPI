% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates Worm BioGRID adjacency matrix.


dt = readtable('../../data/original_data/Worm_BioGRID.txt', 'Delimiter', '\t');
dt = dt( strcmp(table2cell(dt(:,10)), 'taxid:6239')...
    & strcmp(table2cell(dt(:,11)), 'taxid:6239')...
    & strcmp(table2cell(dt(:,12)), 'psi-mi:"MI:0407"(direct interaction)'),:);
id1 = strrep(table2cell(dt(:,1)), 'entrez gene/locuslink:', '');
id2 = strrep(table2cell(dt(:,2)), 'entrez gene/locuslink:', '');

[x, map] = create_matrix(id1, id2, [], 0);
save('../matrix/Worm_BioGRID.mat', 'x', 'map', '-v7.3');