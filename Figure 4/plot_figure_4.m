% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script calculates interaction disorder index (IDI)
% for protein pairs in the golden standards positive (GSP, left panels) and
% and the golden standards negative (GSN, right panels) sets. Computation
% time of CH and AFM


%% Setup the Import Options
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

opts.VariableNames = ["ProteinPairsPos", "AvgINTSPos", "ProteinPairsNeg", "AvgINTSNeg"];
ints = readtable("../INTS_outputs/Yeast_avg_ints_ColabFold_AF2Mv2_unpaired+paired.txt", opts);
opts.VariableNames = ["ProteinPairsPos", "AvgPITMSPos", "ProteinPairsNeg", "AvgPITMSNeg"];
pitms = readtable("../PITMS_outputs/Yeast_avg_pitms_ColabFold_AF2Mv2_unpaired+paired.txt", opts);
ints = table2cell(ints);
pitms = table2cell(pitms);
methods = {'AFM-IS', 'AFM-piTMS', 'CH3-L3'};
% load Yeast DIP network
load('../matrix/Yeast_DIP_net.mat');

% Extract GSP and GSN
limit = 0.2;
num_pairs_GSP_GSN = round(limit * sum(sum(triu(x_lcc,1))));
for i = 1:2
    if i == 1
        P = ints(:,1:2);
        N = ints(:,3:4);
    else
        P = pitms(:,1:2);
        N = pitms(:,3:4);
    end
    
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
    
    scores_AF = [P ; N];
    labels = repelem([1, 0], [length(P) length(N)])';
    % Concatenate the labels vector horizontally to the scores cell
    scores_AF = [scores_AF, num2cell(labels)];
    
    % Sort in descending order the scores cell
    [~,idx] = sortrows(scores_AF, -2);
    sortedScores = scores_AF(idx, :);
    
    TP = {};
    FN = {};
    FP = {};
    TN = {};
    for j = 1:length(sortedScores)
        if sortedScores{j,3} == 1 && j <= length(P)
            TP = [TP ; sortedScores(j, :)];
        elseif sortedScores{j,3} == 1 && j > length(P) 
            FN = [FN ; sortedScores(j, :)];
        elseif sortedScores{j,3} == 0 && j <= length(P)
            FP = [FP ; sortedScores(j, :)];
        elseif sortedScores{j,3} == 0 && j > length(P)
            TN = [TN ; sortedScores(j, :)];
        end
    end
    
    split_array_TP = cellfun(@(x) strsplit(x, '_'), TP(:,1), 'UniformOutput', false);
    split_array_FN = cellfun(@(x) strsplit(x, '_'), FN(:,1), 'UniformOutput', false);
    split_array_FP = cellfun(@(x) strsplit(x, '_'), FP(:,1), 'UniformOutput', false);
    split_array_TN = cellfun(@(x) strsplit(x, '_'), TN(:,1), 'UniformOutput', false);

    %%%%%%%
    C = vertcat(split_array_TP{:});
    mean_avg_iupred_score_TP = zeros(length(C), 1);
    disp([methods{i}, ': TP'])
    for j = 1:length(mean_avg_iupred_score_TP)
        dt_prot_1 = readtable(['../iupred3/results_AF2_Yeast_DIP_Pos/' C{j, 2} '_iupred3.txt']);
        dt_prot_2 = readtable(['../iupred3/results_AF2_Yeast_DIP_Pos/' C{j, 3} '_iupred3.txt']);
        avg_iupred_score_prot_1 = mean(table2array(dt_prot_1(:,3)));
        avg_iupred_score_prot_2 = mean(table2array(dt_prot_2(:,3)));
        
        mean_avg_iupred_score_TP(j) = mean([avg_iupred_score_prot_1, avg_iupred_score_prot_2]);
    end
    
    C = vertcat(split_array_FN{:});
    mean_avg_iupred_score_FN = zeros(length(C), 1);
    disp([methods{i}, ': FN'])
    for j = 1:length(mean_avg_iupred_score_FN)
        dt_prot_1 = readtable(['../iupred3/results_AF2_Yeast_DIP_Pos/' C{j, 2} '_iupred3.txt']);
        dt_prot_2 = readtable(['../iupred3/results_AF2_Yeast_DIP_Pos/' C{j, 3} '_iupred3.txt']);
        avg_iupred_score_prot_1 = mean(table2array(dt_prot_1(:,3)));
        avg_iupred_score_prot_2 = mean(table2array(dt_prot_2(:,3)));
        
        mean_avg_iupred_score_FN(j) = mean([avg_iupred_score_prot_1, avg_iupred_score_prot_2]);
    end
    
    C = vertcat(split_array_FP{:});
    mean_avg_iupred_score_FP = zeros(length(C), 1);
    disp([methods{i}, ': FP'])
    for j = 1:length(mean_avg_iupred_score_FP)
        dt_prot_1 = readtable(['../iupred3/results_AF2_Yeast_DIP_Neg/' C{j, 2} '_iupred3.txt']);
        dt_prot_2 = readtable(['../iupred3/results_AF2_Yeast_DIP_Neg/' C{j, 3} '_iupred3.txt']);
        avg_iupred_score_prot_1 = mean(table2array(dt_prot_1(:,3)));
        avg_iupred_score_prot_2 = mean(table2array(dt_prot_2(:,3)));
        
        mean_avg_iupred_score_FP(j) = mean([avg_iupred_score_prot_1, avg_iupred_score_prot_2]);
    end
    
    C = vertcat(split_array_TN{:});
    mean_avg_iupred_score_TN = zeros(length(C), 1);
    disp([methods{i}, ': TN'])
    for j = 1:length(mean_avg_iupred_score_TN)
        dt_prot_1 = readtable(['../iupred3/results_AF2_Yeast_DIP_Neg/' C{j, 2} '_iupred3.txt']);
        dt_prot_2 = readtable(['../iupred3/results_AF2_Yeast_DIP_Neg/' C{j, 3} '_iupred3.txt']);
        avg_iupred_score_prot_1 = mean(table2array(dt_prot_1(:,3)));
        avg_iupred_score_prot_2 = mean(table2array(dt_prot_2(:,3)));
        
        mean_avg_iupred_score_TN(j) = mean([avg_iupred_score_prot_1, avg_iupred_score_prot_2]);
    end
    
    
    if i == 1
        subplot(4,2,1)
        p = permutation_test(mean_avg_iupred_score_TP, mean_avg_iupred_score_FN, 'median');
        [f,xi,bw1] = ksdensity(mean_avg_iupred_score_TP, 'Support','positive', 'Bandwidth',0.3);
        p1 = plot(xi, f, 'color', 'k');
        % Get the median of the data
        med = median(mean_avg_iupred_score_TP);
        % Add a vertical dashed line at the median
        hold on;
        plot([med, med], [0, max(f)], 'color', 'k');
        hold on;
        [f,xi,bw2] = ksdensity(mean_avg_iupred_score_FN, 'Support','positive', 'Bandwidth',0.3);
        p2 = plot(xi, f, '--', 'color', 'k');
        % Get the median of the data
        med = median(mean_avg_iupred_score_FN);
        % Add a vertical dashed line at the median
        hold on;
        plot([med, med], [0, max(f)], '--', 'color', 'k');
        xlim([0 1]);
        ylabel('Probability');
        if p < 0.001
            title({['GSP = ', num2str(length(P))], 'AFM-IS p-value < 0.001'});
        else
            title({['GSP = ', num2str(length(P))], ['AFM-IS p-value = ' num2str(round(p, 3))]});
        end
        lgd = legend([p1 p2], {'TP','FN'});
        
        subplot(4,2,2)
        p = permutation_test(mean_avg_iupred_score_TN, mean_avg_iupred_score_FP, 'median');
        [f, xi] = ksdensity(mean_avg_iupred_score_TN, 'Support','positive', 'Bandwidth',0.3);
        p1 = plot(xi, f, 'color', 'k');
        % Get the median of the data
        med = median(mean_avg_iupred_score_TN);
        % Add a vertical dashed line at the median
        hold on;
        plot([med, med], [0, max(f)], 'color', 'k');
        hold on;
        [f, xi] = ksdensity(mean_avg_iupred_score_FP, 'Support','positive', 'Bandwidth',0.3);
        p2 = plot(xi, f, '--', 'color', 'k');
        % Get the median of the data
        med = median(mean_avg_iupred_score_FP);
        % Add a vertical dashed line at the median
        hold on;
        plot([med, med], [0, max(f)], '--', 'color', 'k');
        xlim([0 1]);
        if p < 0.001
            title({['GSN = ', num2str(length(N))], 'AFM-IS p-value < 0.001'});
        else
            title({['GSN = ', num2str(length(N))], ['AFM-IS p-value = ' num2str(round(p, 3))]});
        end
        lgd = legend([p1 p2], {'TN','FP'});
        
    elseif i == 2
        subplot(4,2,3)
        p = permutation_test(mean_avg_iupred_score_TP, mean_avg_iupred_score_FN, 'median');
        [f, xi] = ksdensity(mean_avg_iupred_score_TP, 'Support','positive', 'Bandwidth',0.3);
        p1 = plot(xi, f, 'color', 'k');
        % Get the median of the data
        med = median(mean_avg_iupred_score_TP);
        % Add a vertical dashed line at the median
        hold on;
        plot([med, med], [0, max(f)], 'color', 'k');
        hold on;
        [f, xi] = ksdensity(mean_avg_iupred_score_FN, 'Support','positive', 'Bandwidth',0.3);
        p2 = plot(xi, f, '--', 'color', 'k');
        % Get the median of the data
        med = median(mean_avg_iupred_score_FN);
        % Add a vertical dashed line at the median
        hold on;
        plot([med, med], [0, max(f)], '--', 'color', 'k');
        xlim([0 1]);
        ylabel('Probability');
        if p < 0.001
            title('AFM-piTMS p-value < 0.001');
        else
            title(['AFM-piTMS p-value = ' num2str(round(p, 3))]);
        end
        lgd = legend([p1 p2], {'TP','FN'});
        
        subplot(4,2,4)
        p = permutation_test(mean_avg_iupred_score_TN, mean_avg_iupred_score_FP, 'median');
        [f, xi] = ksdensity(mean_avg_iupred_score_TN, 'Support','positive', 'Bandwidth',0.3);
        p1 = plot(xi, f, 'color', 'k');
        % Get the median of the data
        med = median(mean_avg_iupred_score_TN);
        % Add a vertical dashed line at the median
        hold on;
        plot([med, med], [0, max(f)], 'color', 'k');
        hold on;
        [f, xi] = ksdensity(mean_avg_iupred_score_FP, 'Support','positive', 'Bandwidth',0.3);
        p2 = plot(xi, f, '--', 'color', 'k');
        % Get the median of the data
        med = median(mean_avg_iupred_score_FP);
        % Add a vertical dashed line at the median
        hold on;
        plot([med, med], [0, max(f)], '--', 'color', 'k');
        xlim([0 1]);
        if p < 0.001
            title('AFM-piTMS p-value < 0.001');
        else
            title(['AFM-piTMS p-value = ' num2str(round(p, 3))]);
        end
        lgd = legend([p1 p2], {'TN','FP'});
    end    
    
end


load('../network_similarities/CH_L2_L3/results/Yeast_DIP_CHv2_L2_L3');
%load('../network/matrix/Yeast_DIP_net');

CH3_L3_scores = scores_AF;
split_array = cellfun(@(x) strsplit(x, '_'), CH3_L3_scores(:,1), 'UniformOutput', false);
C = vertcat(split_array{:});

S = table2array(scores);
CH_mat = zeros(size(x_lcc, 1), size(x_lcc,2));
ind = sub2ind([size(CH_mat, 1) size(CH_mat, 2)], S(:,1), S(:,2));
CH_mat(ind) = S(:,11); 
CH_mat_sym = CH_mat + CH_mat.';

for i = 1:length(C)
    r_idx = find(ismember(map_lcc, C{i,2}));
    c_idx = find(ismember(map_lcc, C{i,3}));
    CH3_L3_scores{i,2} = CH_mat_sym(r_idx, c_idx);
    
end

% Sort in descending order the scores cell
[~,idx] = sortrows(CH3_L3_scores, -2);
sortedScores = CH3_L3_scores(idx, :);

TP = {};
FN = {};
FP = {};
TN = {};
for i = 1:length(sortedScores)
    if sortedScores{i,3} == 1 && i <= length(P)
        TP = [TP ; sortedScores(i, :)];
    elseif sortedScores{i,3} == 1 && i > length(P)
        FN = [FN ; sortedScores(i, :)];
    elseif sortedScores{i,3} == 0 && i <= length(P)
        FP = [FP ; sortedScores(i, :)];
    elseif sortedScores{i,3} == 0 && i > length(P)
        TN = [TN ; sortedScores(i, :)];
    end
end
        
split_array_TP = cellfun(@(x) strsplit(x, '_'), TP(:,1), 'UniformOutput', false);
split_array_FN = cellfun(@(x) strsplit(x, '_'), FN(:,1), 'UniformOutput', false);
split_array_FP = cellfun(@(x) strsplit(x, '_'), FP(:,1), 'UniformOutput', false);
split_array_TN = cellfun(@(x) strsplit(x, '_'), TN(:,1), 'UniformOutput', false);

C = vertcat(split_array_TP{:});
mean_avg_iupred_score_TP = zeros(length(C), 1);
disp('CH3-L3 TP')
for i = 1:length(mean_avg_iupred_score_TP)
    dt_prot_1 = readtable(['../iupred3/results_AF2_Yeast_DIP_Pos/' C{i, 2} '_iupred3.txt']);
    dt_prot_2 = readtable(['../iupred3/results_AF2_Yeast_DIP_Pos/' C{i, 3} '_iupred3.txt']);
    avg_iupred_score_prot_1 = mean(table2array(dt_prot_1(:,3)));
    avg_iupred_score_prot_2 = mean(table2array(dt_prot_2(:,3)));
    
    mean_avg_iupred_score_TP(i) = mean([avg_iupred_score_prot_1, avg_iupred_score_prot_2]);
end

C = vertcat(split_array_FN{:});
mean_avg_iupred_score_FN = zeros(length(C), 1);
disp('CH3-L3 FN')
for i = 1:length(mean_avg_iupred_score_FN)
    dt_prot_1 = readtable(['../iupred3/results_AF2_Yeast_DIP_Pos/' C{i, 2} '_iupred3.txt']);
    dt_prot_2 = readtable(['../iupred3/results_AF2_Yeast_DIP_Pos/' C{i, 3} '_iupred3.txt']);
    avg_iupred_score_prot_1 = mean(table2array(dt_prot_1(:,3)));
    avg_iupred_score_prot_2 = mean(table2array(dt_prot_2(:,3)));
    
    mean_avg_iupred_score_FN(i) = mean([avg_iupred_score_prot_1, avg_iupred_score_prot_2]);
end

C = vertcat(split_array_FP{:});
mean_avg_iupred_score_FP = zeros(length(C), 1);
disp('CH3-L3 FP')
for i = 1:length(mean_avg_iupred_score_FP)
    dt_prot_1 = readtable(['../iupred3/results_AF2_Yeast_DIP_Neg/' C{i, 2} '_iupred3.txt']);
    dt_prot_2 = readtable(['../iupred3/results_AF2_Yeast_DIP_Neg/' C{i, 3} '_iupred3.txt']);
    avg_iupred_score_prot_1 = mean(table2array(dt_prot_1(:,3)));
    avg_iupred_score_prot_2 = mean(table2array(dt_prot_2(:,3)));
    
    mean_avg_iupred_score_FP(i) = mean([avg_iupred_score_prot_1, avg_iupred_score_prot_2]);
end

C = vertcat(split_array_TN{:});
mean_avg_iupred_score_TN = zeros(length(C), 1);
disp('CH3-L3 TN')
for i = 1:length(mean_avg_iupred_score_TN)
    dt_prot_1 = readtable(['../iupred3/results_AF2_Yeast_DIP_Neg/' C{i, 2} '_iupred3.txt']);
    dt_prot_2 = readtable(['../iupred3/results_AF2_Yeast_DIP_Neg/' C{i, 3} '_iupred3.txt']);
    avg_iupred_score_prot_1 = mean(table2array(dt_prot_1(:,3)));
    avg_iupred_score_prot_2 = mean(table2array(dt_prot_2(:,3)));
    
    mean_avg_iupred_score_TN(i) = mean([avg_iupred_score_prot_1, avg_iupred_score_prot_2]);
end


subplot(4,2,5)
p = permutation_test(mean_avg_iupred_score_TP, mean_avg_iupred_score_FN, 'median');
[f, xi] = ksdensity(mean_avg_iupred_score_TP, 'Support','positive', 'Bandwidth',0.3);
p1 = plot(xi, f, 'color', 'r');
% Get the median of the data
med = median(mean_avg_iupred_score_TP);
% Add a vertical dashed line at the median
hold on;
plot([med, med], [0, max(f)], 'color', 'r');
hold on;
[f, xi] = ksdensity(mean_avg_iupred_score_FN, 'Support','positive', 'Bandwidth',0.3);
p2 = plot(xi, f, '--', 'color', 'r');
% Get the median of the data
med = median(mean_avg_iupred_score_FN);
% Add a vertical dashed line at the median
hold on;
plot([med, med], [0, max(f)], '--', 'color', 'r');
xlim([0 1]);
xlabel('Interaction Disorder Index');
ylabel('Probability');        
if p < 0.001
    title('CH3-L3 p-value < 0.001');
else
    title(['CH3-L3 p-value = ' num2str(round(p, 3))]);
end
lgd = legend([p1 p2], {'TP','FN'});

subplot(4,2,6)
p = permutation_test(mean_avg_iupred_score_TN, mean_avg_iupred_score_FP, 'median');
[f, xi] = ksdensity(mean_avg_iupred_score_TN, 'Support','positive', 'Bandwidth',0.3);
p1 = plot(xi, f, 'color', 'r');
% Get the median of the data
med = median(mean_avg_iupred_score_TN);
% Add a vertical dashed line at the median
hold on;
plot([med, med], [0, max(f)], 'color', 'r');
hold on;
[f, xi] = ksdensity(mean_avg_iupred_score_FP, 'Support','positive', 'Bandwidth',0.3);
p2 = plot(xi, f, '--', 'color', 'r');
% Get the median of the data
med = median(mean_avg_iupred_score_FP);
% Add a vertical dashed line at the median
hold on;
plot([med, med], [0, max(f)], '--', 'color', 'r');
xlim([0 1]);
xlabel('Interaction Disorder Index');
if p < 0.001
    title('CH3-L3 p-value < 0.001');
else
    title(['CH3-L3 p-value = ' num2str(round(p, 3))]);
end
lgd = legend([p1 p2], {'TN','FP'});


%%%%%% CH and AFM time calculation %%%%%%
disp('CH and AFM time calculation')
output_dir_pos = '../../data/Yeast_AFM_output/Positive_set/AlphaFold2-multimer-v2/unpaired+paired/';
output_dir_neg = '../../data/Yeast_AFM_output/Negative_set/AlphaFold2-multimer-v2/unpaired+paired/';

load('../matrix/list_pairs_1percent_Yeast.mat');

methods = {'RA_L2','CH1_L2','CH2_L2','CH3_L2','RA_L3','CH1_L3','CH2_L3','CH3_L3'};
CHA_option = {'CH2_L2', 'CH3_L2', 'CH2_L3', 'CH3_L3'};

% Create a logical mask for the upper triangle of x
mask = triu(ones(size(x_lcc)), 1);

% Find the indices of missing links (where x is 0) in the upper triangle
[row, col] = find(x_lcc == 0 & mask);

% Create the w matrix with missing link indices
w = [row, col];

% Run CHA and save elapsed time. Run with parallelization on 12 workers.
tic;
%[scores, CHA_info] = CHA_linkpred_monopartite(x_lcc, methods, CHA_option);
CH_Ln_mex(x_lcc, 3, w, 1);
CH_total_elapsed_time = toc;

% Calculate number of missing links in Yeast net and CH average time for
% one link.
num_missing_links = sum(triu(x_lcc == 0, 1), 'all');
CH_time_one_link = CH_total_elapsed_time / num_missing_links;
% 0.013ms

% Concatenate positive and negative sets and extract unique pairs.
concatenate_L_pairs = L_pairs(:);
unique_pairs = unique(vertcat(concatenate_L_pairs{:}), 'rows');

% Convert node indices in the net to their respective uniprot ID.
uniprot_ids = map_lcc(unique_pairs);

% Get list of folders 
struct_folders_neg = dir([output_dir_neg, '*.fasta']);
struct_folders_pos = dir([output_dir_pos, '*.fasta']);
list_folders = {struct_folders_neg.name, struct_folders_pos.name}';

fileName = 'timings.json';
total_time_in_s = zeros(length(uniprot_ids), 1);
count = 0;
for i = 1:length(uniprot_ids)
    idx = find(contains(list_folders, [uniprot_ids{i,1}, '_', uniprot_ids{i,2}]));
    if isempty(idx)
        idx = find(contains(list_folders, [uniprot_ids{i,2}, '_', uniprot_ids{i,1}]));
        if isempty(idx)
            %disp(i)
            count = count + 1;
            total_time_in_s(i) = 0;
            continue;
        end
    end
    if idx <= length(struct_folders_neg)
        path_timings = [output_dir_neg, list_folders{idx}, '\' ,fileName];
    else
        path_timings = [output_dir_pos, list_folders{idx}, '\' ,fileName];
        
    end
    fid = fopen(path_timings); % Opening the file
    raw = fread(fid,inf); % Reading the contents
    str = char(raw'); % Transformation
    fclose(fid); % Closing the file
    data = struct2cell(jsondecode(str)); % Using the jsondecode function to parse JSON from string
    time_in_s = sum([data{:}]);
    total_time_in_s(i) = time_in_s;
    
end

subplot(4,2,7);
AFM_time_one_link = sum(total_time_in_s)/(length(uniprot_ids) - count);
b = bar(log10([CH_time_one_link AFM_time_one_link]), 'FaceColor','flat');
b.CData(1,:) = [1 0 0];
b.CData(2,:) = [0 0 0];
ylim([-6 6]);
set(gca,'xticklabel',{'CH' 'AFM'});
ylabel('log10(time in s)');
ctr(1,:) = bsxfun(@plus, b(1).XData, b(1).XOffset');   
ydt(1,:) = b(1).YData;  
CH_one_link_relevant_time = round(CH_time_one_link * 10^6);
AFM_one_link_relevant_time = round(AFM_time_one_link / 60);
text(ctr(1,:), ydt(1,:) - 0.06, {[num2str(CH_one_link_relevant_time), ' µs']; [num2str(AFM_one_link_relevant_time), ' min'] }, 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',12, 'Color','k')

subplot(4,2,8);
AFM_total_time_missing_links = AFM_time_one_link * num_missing_links;
b = bar(log10([CH_total_elapsed_time AFM_total_time_missing_links]), 'FaceColor', 'flat');
b.CData(1,:) = [1 0 0];
b.CData(2,:) = [0 0 0];
ylim([0 12]);
set(gca,'xticklabel',{'CH' 'AFM'});
ctr(1,:) = bsxfun(@plus, b(1).XData, b(1).XOffset');   
ydt(1,:) = b(1).YData;  
% In one century, there are 3.1536 giga seconds approximately
AFM_total_relevant_time_missing_links = round(AFM_total_time_missing_links / 3.1536e9, 1);
text(ctr(1,:), ydt(1,:) - 0.06, {[num2str(round(CH_total_elapsed_time / 60, 1)), ' min']; [num2str(AFM_total_relevant_time_missing_links), ' cent'] }, 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',12, 'Color','k')

