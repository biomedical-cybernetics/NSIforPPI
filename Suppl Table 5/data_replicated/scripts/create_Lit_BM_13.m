% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates Lit-BM-13 adjacency matrix.


dt = readtable('../../data/original_data/Lit-BM-13.tsv', 'FileType', 'text','Delimiter', '\t');
dt = table2array(dt(:, [1 3]));

[x, map] = create_matrix(dt(:,1), dt(:,2), [], 0);
save('../matrix/Lit_BM_13.mat', 'x', 'map', '-v7.3');
