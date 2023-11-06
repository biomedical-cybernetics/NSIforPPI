% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script conducts 10 simulations, each time
% generating a list of positive and negative protein interactions.
% These interactions are randomly selected from the positive and negative
% sets, respectively, and their quantity is set to 1% of the total links in
% the network. The script then removes the positively interacting protein
% pairs from the network, causing a perturbation. As output, the script
% returns 10 lists and 10 perturbed networks.


path_network = '../matrix/';
load([path_network, 'Yeast_DIP_net.mat']);
pos_folder = dir('../../data/FASTA_files/Yeast/Yeast_Positive_set/*.fasta');
pos_filenames = {pos_folder.name}';
P = regexp(pos_filenames, '[_.]', 'split');
P = vertcat(P{:});
P = P(:,2:3);

neg_folder = dir('../../data/FASTA_files/Yeast/Yeast_Negative_set/*.fasta');
neg_filenames = {neg_folder.name}';
N = regexp(neg_filenames, '[_.]', 'split');
N = vertcat(N{:});
N = N(:,2:3);

idx1_P = zeros(length(P), 1);
idx2_P = zeros(length(P), 1);
idx1_N = zeros(length(N), 1);
idx2_N = zeros(length(N), 1);
%assert(length(P) == length(N), 'Positive and negative sets have different size.');

for k = 1:length(P)
    idx1_P(k) = find(strcmp(map_lcc, P{k,1}));
    idx2_P(k) = find(strcmp(map_lcc, P{k,2}));
    
end

for k = 1:length(N)
    idx1_N(k) = find(strcmp(map_lcc, N{k,1}));
    idx2_N(k) = find(strcmp(map_lcc, N{k,2}));
    
end

full_set = cat(1,[idx1_P, idx2_P], [idx1_N, idx2_N]);
n = 0;
r = 1;
list_seeds = zeros(10, 1);
while r<=10
    rng(n);
    x_rem = random_link_removal(x_lcc, 0.01, 'rand_conn');
    [n1,n2] = find(triu(x_lcc-x_rem,1));
    count = 0;
    for i = 1:length(n1)
        if ~(ismember([n1(i), n2(i)], full_set, 'row') || ismember([n2(i), n1(i)], full_set, 'row'))
            %disp(i);
            count = count + 1;
        end
    end
    if count == 0
        disp(n);
        list_seeds(r) = n;
        r = r + 1;
        n = n + 1;
    else
        n = n + 1;
    end
end

% positive pairs in first column
% negative pairs in second column
L_pairs = cell(10, 2);
L_net = cell(10, 1);
for i = 1:length(list_seeds)
    rng(list_seeds(i));
    x_rem = random_link_removal(x_lcc, 0.01, 'rand_conn');
    [n1,n2] = find(triu(x_lcc-x_rem,1));
    
    rng(list_seeds(i));
    idx = randi([1 length(N)],1, length(n1))';
    L_pairs(i,:) = {[n1, n2], [idx1_N(idx), idx2_N(idx)]};
    L_net{i} = x_rem;
    
end

save('../matrix/network_perturbed_1percent_Yeast.mat', 'L_net');
save('../matrix/list_pairs_1percent_Yeast.mat', 'L_pairs');
