% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates the GSP and GSN sets and perturbed
% the network by deleting 5% of the links.


opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["ProteinPairsPos", "AvgPairwiseRMSDPos", "ProteinPairsNeg", "AvgPairwiseRMSDNeg"];
opts.VariableTypes = ["char", "double", "char", "double"];
opts = setvaropts(opts, [1, 3], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 3], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
opts.VariableNames = ["ProteinPairsPos", "AvgINTSPos", "ProteinPairsNeg", "AvgINTSNeg"];
ints = readtable("../INTS_outputs/Yeast_avg_ints_ColabFold_AF2Mv2_unpaired+paired.txt", opts);
opts.VariableNames = ["ProteinPairsPos", "AvgPITMSPos", "ProteinPairsNeg", "AvgPITMSNeg"];
pitms = readtable("../PITMS_outputs/Yeast_avg_pitms_ColabFold_AF2Mv2_unpaired+paired.txt", opts);
load('../matrix/Yeast_DIP_net.mat');

%% Convert to output type
ints = table2cell(ints);
pitms = table2cell(pitms);
%numIdx = cellfun(@(x) ~isnan(str2double(x)), ints);
%ints(numIdx) = cellfun(@(x) {str2double(x)}, ints(numIdx));
%numIdx = cellfun(@(x) ~isnan(str2double(x)), pitms);
%pitms(numIdx) = cellfun(@(x) {str2double(x)}, pitms(numIdx));


%% Clear temporary variables
clear opts

% Top 20% positive set
% Bottom 20% negative set
limit = 0.2;
num_pairs_GSP_GSN = round(limit * sum(sum(triu(x_lcc,1))));

P = ints(:,1:2);
N = ints(:,3:4);
% Identify the rows with NaN values
nanRowsP = any(cellfun(@(x) any(isnan(x)), P), 2);
nanRowsN = any(cellfun(@(x) any(isnan(x)), N), 2);

% Remove the rows with NaN values
P(nanRowsP, :) = [];
N(nanRowsN, :) = [];

% Extract the numeric portion from the first column
splitValues_P = cellfun(@(x) strsplit(x, '_'), P(:, 1), 'UniformOutput', false);
splitValues_N = cellfun(@(x) strsplit(x, '_'), N(:, 1), 'UniformOutput', false);
numbers_P = cellfun(@(x) str2double(x(1)), splitValues_P);
numbers_N = cellfun(@(x) str2double(x(1)), splitValues_N);
%numbers_2 = cellfun(@(x) x(1), splitValues);

% Sort the numbers in ascending order and get the sorting indices
[sortedNumbers, sortingIndices_P] = sort(numbers_P);
[sortedNumbers, sortingIndices_N] = sort(numbers_N);

% Sort the cell array based on the sorting indices and select GSP and
% GSN.
P = P(sortingIndices_P(1:num_pairs_GSP_GSN), :);
N = N(sortingIndices_N(1:num_pairs_GSP_GSN), :);

tmp_P = regexp(P(:,1), '[_.]', 'split');
tmp_N = regexp(N(:,1), '[_.]', 'split');
tmp_P = vertcat(tmp_P{:});
tmp_N = vertcat(tmp_N{:});
P = tmp_P(:,2:3);
N = tmp_N(:,2:3);

idx1_P = zeros(length(P), 1);
idx2_P = zeros(length(P), 1);
idx1_N = zeros(length(N), 1);
idx2_N = zeros(length(N), 1);

for k = 1:length(P)
    node1 = find(strcmp(map_lcc, P{k,1}));
    node2 = find(strcmp(map_lcc, P{k,2}));
    % iterative exclusion of leaves
    if sum(x_lcc(node1, :)) > 1 && sum(x_lcc(node2, :)) > 1
        idx1_P(k) = node1;
        idx2_P(k) = node2;
    end
    
end
idx1_P = idx1_P(idx1_P > 0); 
idx2_P = idx2_P(idx2_P > 0); 

for k = 1:length(N)
    idx1_N(k) = find(strcmp(map_lcc, N{k,1}));
    idx2_N(k) = find(strcmp(map_lcc, N{k,2}));
    
end

s = RandStream("dsfmt19937");
n_sims = 10;
n_percent = 5;
n_pairs_to_remove = round(sum(triu(x_lcc, 1), 'all') * (n_percent/100)); % number of links to delete
%idx_pos_pairs = randi(s, [1 length(P)], n_sims, n_pairs_to_remove)';
%idx_neg_pairs = randi(s, [1 length(N)], n_sims, n_pairs_to_remove)';

% take upper triangle
x = triu(x_lcc,1);

d = sum(sum(x));            % number of existing links

% positive pairs in first column
% negative pairs in second column
L_pairs = cell(10, 2);
L_net = cell(10, 1);
ncc = Inf;
iter = 0;
while iter < n_sims
    temp = x;
    
    % randomly sample links to delete
    %idx = randsample(dr, num);
    i_pos = randsample(s, length(idx1_P), n_pairs_to_remove);
    i_neg = randsample(s, length(N), n_pairs_to_remove);
    
    % delete links
    temp(sub2ind(size(temp), idx1_P(i_pos,1), idx2_P(i_pos,1))) = 0;
    temp(sub2ind(size(temp), idx2_P(i_pos,1), idx1_P(i_pos,1))) = 0;
    %temp(sub2ind(size(temp), idx1_N(i_neg,1), idx2_N(i_neg,1))) = 0;
    
    % compute the number of connected components
    % for undirected graphs the algorithm considers only the lower triangle
    [ncc,~] = graphconncomp(temp', 'Directed', false);
    %disp(ncc);
    if ncc < 2
        iter = iter + 1;
        disp(iter);
        % return the matrix symmetric
        %x = temp; clear temp
        [n1,n2] = find(triu(x-temp,1));
        L_net{iter} = max(temp, temp');
        L_pairs(iter,:) = {[n1, n2], [idx1_N(i_neg, 1), idx2_N(i_neg, 1)]};
    end
    
end


save('../matrix/network_perturbed_5percent_GSP_GSN_Yeast.mat', 'L_net');
save('../matrix/list_pairs_5percent_GSP_GSN_Yeast.mat', 'L_pairs');
