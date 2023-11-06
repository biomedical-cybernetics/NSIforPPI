% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates Arabidopsis adjacency matrix.


dt = readtable('../../data/original_data/Arabidopsis.xls');
dt = table2array(dt(:, [1 2]));

[x, map] = create_matrix(dt(:,1), dt(:,2), [], 0);
save('../matrix/Arabidopsis.mat', 'x', 'map', '-v7.3');