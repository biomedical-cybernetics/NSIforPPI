% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates Bioplex adjacency matrix.


dt = readtable('../../data/original_data/Bioplex.tsv', 'FileType', 'text','Delimiter', '\t');
dt = table2array(dt(:, [1 2]));

[x, map] = create_matrix(dt(:,1), dt(:,2), [], 0);
save('../matrix/Bioplex.mat', 'x', 'map', '-v7.3');