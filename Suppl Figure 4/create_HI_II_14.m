% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates HI-II-14 adjacency matrix.


dt = readtable('../../data/original_data/HI_II_14.txt', 'Delimiter', '\t');
dt = dt( ~strcmp(table2cell(dt(:,1)), '-')...
    & ~strcmp(table2cell(dt(:,2)), '-'),1:2);
id1 = strrep(table2cell(dt(:,1)), 'uniprotkb:', '');
id2 = strrep(table2cell(dt(:,2)), 'uniprotkb:', '');

[x, map] = create_matrix(id1, id2, [], 0);
save('../matrix/HI_II_14.mat', 'x', 'map', '-v7.3');
