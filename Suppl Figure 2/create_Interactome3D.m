% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates Interactome3D adjacency matrix.


adj_list = readtable('../../data/original_data/interactions2015_12.dat', 'Delimiter', '\t');
adj_list = table2array(adj_list(:, [1 2]));

[x, map] = create_matrix(adj_list(:,1), adj_list(:,2), [], 0);
save('../matrix/Interactome3D.mat', 'x', 'map', '-v7.3');
