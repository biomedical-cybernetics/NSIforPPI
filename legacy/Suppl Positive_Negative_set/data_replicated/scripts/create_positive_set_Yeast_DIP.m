% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2022
% Description: This MATLAB script creates the positive set of the Yeast DIP
% network. The positive set is composed by all the existing interactions 
% with a Wang GO semantic similarity. 


file_path_net = '../matrix/';
file_path_generated_data = '../';
organism = 'Yeast';
load([file_path_net organism '_DIP_net.mat']);

% read conversion file
data_conversion = readcell([file_path_generated_data 'Uniprot_to_EntrezID_' organism '_DIP_net.txt'], 'NumHeaderLines', 1);
M_BP = readmatrix([file_path_generated_data organism '_DIP_BP.txt']);
% Extract Entrez ID in BP Wang GO semantic similarity matrix
ids_BP = M_BP(:,1);
M_BP(:,1) = [];
assert(size(M_BP, 1) == size(M_BP, 2), 'Matrix is not squared.')
% Extract Entrez ID in CC Wang GO semantic similarity matrix
M_CC = readmatrix([file_path_generated_data organism '_DIP_CC.txt']);
ids_CC = M_CC(:,1);
M_CC(:,1) = [];
assert(size(M_CC, 1) == size(M_CC, 2), 'Matrix is not squared.')


% find elements of 'ids_CC' in 'ids_BP' and vice versa and return
% corresponding indices.
ids_to_keep_CC = find(ismember(ids_CC, ids_BP));
ids_to_keep_BP = find(ismember(ids_BP, ids_CC));
% Ensure 'M_CC' and 'M_BP' have the same size.
M_CC = M_CC(ids_to_keep_CC, ids_to_keep_CC);
M_BP = M_BP(ids_to_keep_BP, ids_to_keep_BP);
assert(any([size(M_CC,1) size(M_CC,2)] == [size(M_BP,1) size(M_BP,2)]), 'M_BP and M_CC not same dimension.');
% Take min values between 'M_BP' and 'M_CC'.
res_M = min(M_BP, M_CC);
assert(isequal(ids_BP(ids_to_keep_BP), ids_CC(ids_to_keep_CC)), 'ids rowwise and columnwise not identical.')
ids_res_M = ids_BP(ids_to_keep_BP);
% Sort in descending order values in upper triangle of 'res_M'.
% Replace the lower triangle, including the diagonal, with -Inf 
% to exclude these elements from the sorting operation.
UI = tril(true(size(res_M)));
res_M(UI) = -Inf;
% Sort in descending order values in upper triangle of 'res_M'.
[out, indices] = sort(res_M(:),'descend');
% remove -Inf values.
inf_idx = find(out == -Inf);
indices(inf_idx) = [];
% find indices of sorted elements.
[row, col] = ind2sub([size(res_M, 1) size(res_M, 2)], indices);

% An Entrez ID may have multiple Uniprot IDs
[unique_entrez_ids, ia] = unique([data_conversion{:,2}]','stable');
data_conversion = data_conversion(ia,:);
% while selecting the top 5% pairs in the distribution of 'res_M',
% they should be interacting pairs in the network 'x_lcc'.

% Determine the count of interacting protein pairs within the Yeast network
% for which a Wang Gene Ontology (GO) semantic similarity has been computed
ii = find(ismember([data_conversion{:,2}]', ids_res_M));
ii = find(ismember(map_lcc, data_conversion(ii,1)));
if length(ii) ~= length(ids_res_M)
    error('One or multiple Uniprot ID(s) of map_lcc not found in data_conversion');
end
total_interacting_pairs = sum(sum(triu(x_lcc(ii,ii))));

interacting_pairs = 0;
pairs_to_test = {'Protein 1' 'Protein 2'};
num_unique_pairs = sum(triu(true(size(res_M)),1), 'all');
for i = 1:round(num_unique_pairs)
    %disp(i);
    idx_ids_r = find([data_conversion{:,2}]' == ids_res_M(row(i)));
    row_id_x = find(ismember(map_lcc, data_conversion{idx_ids_r(1),1}));
    idx_ids_c = find([data_conversion{:,2}]' == ids_res_M(col(i)));
    col_id_x = find(ismember(map_lcc, data_conversion{idx_ids_c(1),1}));
    
    if x_lcc(row_id_x, col_id_x) == 1
        pairs_to_test = [pairs_to_test; [map_lcc(row_id_x) map_lcc(col_id_x)]];
        interacting_pairs = interacting_pairs + 1;
        disp(['interacting pairs : ', num2str(interacting_pairs)]);
        if interacting_pairs == total_interacting_pairs
            break;
        end
    end
    
end

writecell(pairs_to_test, [file_path_generated_data 'Positive_set_' organism '.txt'], 'Delimiter','tab');
