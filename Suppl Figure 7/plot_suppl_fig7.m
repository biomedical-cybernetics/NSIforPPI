% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script generates the performance of the AFM settings in predicting
% the GSP and GSN protein pairs of the Yeast DIP network from original data.





%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: Z:\Master students\Ilyes Abdelhamid\AF2\data\alphafold_output\Yeast_avg_pairwise_rmsd_AF2_align.txt
%
% Auto-generated by MATLAB on 11-Feb-2022 09:18:20
 
%% Setup the Import Options
tic;
AF2_variants = {'AF2Mv1_unpaired+paired', 'AF2Mv1_unpaired', 'AF2Mv1_paired', ...
    'AF2Mv2_unpaired+paired', 'AF2Mv2_unpaired', 'AF2Mv2_paired'};
methods = {'interface score', 'piTM score'};
    res_figure_AUC_mROC = cell(length(AF2_variants), length(methods));
    res_figure_AUC_ROC = cell(length(AF2_variants), length(methods));
    res_figure_AUC_PR = cell(length(AF2_variants), length(methods));
    res_figure_NDCG = cell(length(AF2_variants), length(methods));
    res_figure_MCC = cell(length(AF2_variants), length(methods));
for i = 1:length(AF2_variants)
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
    ints = readtable(['../INTS_outputs/Yeast_avg_ints_ColabFold_', AF2_variants{i}, '.txt'], opts);
    opts.VariableNames = ["ProteinPairsPos", "AvgPITMSPos", "ProteinPairsNeg", "AvgPITMSNeg"];
    pitms = readtable(['../PITMS_outputs/Yeast_avg_pitms_ColabFold_', AF2_variants{i}, '.txt'], opts);
    
    load('../matrix/Yeast_DIP_net.mat');
    
    %% Convert to output type
    ints = table2cell(ints);
    pitms = table2cell(pitms);
    numIdx = cellfun(@(x) ~isnan(str2double(x)), ints);
    ints(numIdx) = cellfun(@(x) {str2double(x)}, ints(numIdx));
    numIdx = cellfun(@(x) ~isnan(str2double(x)), pitms);
    pitms(numIdx) = cellfun(@(x) {str2double(x)}, pitms(numIdx));
    
    %% Clear temporary variables
    clear opts
    disp(AF2_variants{i});
    steps = linspace(0.01, 0.2, 20);
    
    tmp_res_figure_AUC_mROC = zeros(20, 2);
    tmp_res_figure_AUC_ROC = zeros(20, 2);
    tmp_res_figure_AUC_PR = zeros(20, 2);
    tmp_res_figure_NDCG = zeros(20, 2);
    tmp_res_figure_MCC = zeros(20, 2);
    for j = 1:length(steps)
        threshold = steps(j);
        num_pairs_to_test = round(sum(triu(x_lcc,1),'all') * threshold);
        
        
        C_ints_pos = regexp(ints(:,1), '_', 'split');
        C_ints_neg = regexp(ints(:,3), '_', 'split');
        index_pos = cellfun(@(x) ~all(strcmp(x, 'None')), C_ints_pos, 'UniformOutput', 1);
        index_neg = cellfun(@(x) ~all(strcmp(x, 'None')), C_ints_neg, 'UniformOutput', 1);
        C_ints_pos = vertcat(C_ints_pos{index_pos});
        C_ints_neg = vertcat(C_ints_neg{index_neg});
        column1 = str2double(C_ints_pos(:,1));
        column3 = str2double(C_ints_neg(:,1));
        [~, sortOrder_pos] = sort(column1);
        [~, sortOrder_neg] = sort(column3);
        ints_pos = ints(sortOrder_pos(1:num_pairs_to_test), 1:2);
        ints_neg = ints(sortOrder_neg(1:num_pairs_to_test), 3:4);
        
        C_pitms_pos = regexp(pitms(:,1), '_', 'split');
        C_pitms_neg = regexp(pitms(:,3), '_', 'split');
        index_pos = cellfun(@(x) ~all(strcmp(x, 'None')), C_pitms_pos, 'UniformOutput', 1);
        index_neg = cellfun(@(x) ~all(strcmp(x, 'None')), C_pitms_neg, 'UniformOutput', 1);
        C_pitms_pos = vertcat(C_pitms_pos{index_pos});
        C_pitms_neg = vertcat(C_pitms_neg{index_neg});
        column1 = str2double(C_pitms_pos(:,1));
        column3 = str2double(C_pitms_neg(:,1));
        [~, sortOrder_pos] = sort(column1);
        [~, sortOrder_neg] = sort(column3);
        pitms_pos = pitms(sortOrder_pos(1:num_pairs_to_test), 1:2);
        pitms_neg = pitms(sortOrder_neg(1:num_pairs_to_test), 3:4);
        
        labels = repelem([1, 0], num_pairs_to_test)';
        predicted_scores_INTS = [ints_pos(:,2); ints_neg(:,2)];
        predicted_scores_PITMS = [pitms_pos(:,2); pitms_neg(:,2)];
        
        [measures_INTS, ~, ~] = prediction_evaluation(cell2mat(predicted_scores_INTS), labels);
        [measures_PITMS, ~, ~] = prediction_evaluation(cell2mat(predicted_scores_PITMS), labels);
        tmp_res_figure_AUC_mROC(j,:) = [measures_INTS.auc_mroc ; measures_PITMS.auc_mroc];
        tmp_res_figure_AUC_ROC(j,:) = [measures_INTS.auc_roc ; measures_PITMS.auc_roc];
        tmp_res_figure_AUC_PR(j,:) = [measures_INTS.auc_pr ; measures_PITMS.auc_pr];
        tmp_res_figure_NDCG(j,:) = [measures_INTS.ndcg ; measures_PITMS.ndcg];
        tmp_res_figure_MCC(j,:) = [measures_INTS.mcc ; measures_PITMS.mcc];
        
    end
    res_figure_AUC_mROC(i, :) = {tmp_res_figure_AUC_mROC(:,1), tmp_res_figure_AUC_mROC(:,2)}; 
    res_figure_AUC_ROC(i, :) = {tmp_res_figure_AUC_ROC(:,1), tmp_res_figure_AUC_ROC(:,2)}; 
    res_figure_AUC_PR(i, :) = {tmp_res_figure_AUC_PR(:,1), tmp_res_figure_AUC_PR(:,2)};
    res_figure_NDCG(i, :) = {tmp_res_figure_NDCG(:,1), tmp_res_figure_NDCG(:,2)}; 
    res_figure_MCC(i, :) = {tmp_res_figure_MCC(:,1), tmp_res_figure_MCC(:,2)}; 
    
end

%%% AUC-ROC BAR PLOT %%%
subplot(2,1,1);
Y = cell2mat(cellfun(@mean,res_figure_AUC_ROC,'uni',0));
stderror_series = cell2mat(cellfun(@std,res_figure_AUC_ROC,'uni',0))/sqrt(size(res_figure_AUC_ROC, 1));

b = bar(Y);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.7 0.7 0.7];
                                                 % Return �bar� Handle
for k1 = 1:size(res_figure_AUC_ROC, 2)
    ctr(k1,:) = bsxfun(@plus, b(1).XData, b(k1).XOffset');    % Note: �XOffset� Is An Undocumented Feature, This Selects The �bar� Centres
    ydt(k1,:) = b(k1).YData;                                     % Individual Bar Heights
    text(ctr(k1,:), ydt(k1,:) + 0.02, sprintfc('%.2f', ydt(k1,:)), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',12, 'Color','k')
end

set(gca, 'XTickLabel', '');
ylim([0.5 1])
ylabel('AUC-ROC')
box off

hold on
err = errorbar(ctr, ydt, stderror_series', 'r', 'linestyle', 'none');
for i = 1:length(AF2_variants)
    err(i).LineWidth = 1.5;
end
box off
hold off

%%% AUC-PR BAR PLOT %%%
subplot(2,1,2);
Y = cell2mat(cellfun(@mean,res_figure_AUC_PR,'uni',0));
stderror_series = cell2mat(cellfun(@std,res_figure_AUC_PR,'uni',0))/sqrt(size(res_figure_AUC_PR, 1));

b = bar(Y);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.7 0.7 0.7];
                                                 
for k1 = 1:size(res_figure_AUC_PR, 2)
    ctr(k1,:) = bsxfun(@plus, b(1).XData, b(k1).XOffset');    
    ydt(k1,:) = b(k1).YData;                                  
    text(ctr(k1,:), ydt(k1,:) + 0.02, sprintfc('%.2f', ydt(k1,:)), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',12, 'Color','k')
end

xl = {'unpaired+paired (AFM-v1)', 'unpaired (AFM-v1)', ...
    'paired (AFM-v1)', 'unpaired+paired (AFM-v2)', ...
    'unpaired (AFM-v2)', 'paired (AFM-v2)'};
xl = cellfun(@(x) strrep(x,' ','\newline'), xl,'UniformOutput',false);
set(gca, 'XTickLabel', xl, 'FontSize', 12);
ylim([0.5 1])
ylabel('AUC-PR')
box off

hold on
err = errorbar(ctr, ydt, stderror_series', 'r', 'linestyle', 'none');
for i = 1:length(AF2_variants)
    err(i).LineWidth = 1.5;
end
box off
hold off
toc;
