% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates Yeast adjacency matrix.


dt = readtable('../../data/original_data/Yeast.txt', 'Delimiter', '\t', 'ReadVariableName', false);
dt = table2array(dt(:, [1 2]));
[x, map] = create_matrix(dt(:,1), dt(:,2), [], 0);

save('../matrix/Yeast.mat', 'x', 'map', '-v7.3');