% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates Hein adjacency matrix.


dt = readtable('../../data/original_data/Hein.xlsx');
dt = table2array(dt(:, 3:4));
id1 = cell(size(dt,1),1);
id2 = cell(size(dt,1),1);
for i = 1:size(dt,1)
    if length(strsplit(char(dt{i,1}), ';')) > 1
        tmp_string = strsplit(char(dt{i,1}), ';');
        id1{i} = tmp_string{1};
    else
        id1{i} = table2array(dt(i,1));
    end
    if length(strsplit(char(dt{i,2}), ';')) > 1
        tmp_string = strsplit(char(dt{i,2}), ';');
        id2{i} = tmp_string{1};
    else
        id2{i} = table2array(dt(i,2));
    end
end

adj_list = [id1 id2];

[x, map] = create_matrix(id1, id2, [], 0);
save('../matrix/Hein.mat', 'x', 'map', '-v7.3');









