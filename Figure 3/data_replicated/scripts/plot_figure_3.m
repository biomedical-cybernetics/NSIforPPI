% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script plots Figure 3 with the replicated results.


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
 
% Import the data
ints = readtable("../INTS_outputs/Yeast_avg_ints_ColabFold_AF2Mv2_unpaired+paired.txt", opts);
opts.VariableNames = ["ProteinPairsPos", "AvgPITMSPos", "ProteinPairsNeg", "AvgPITMSNeg"];
pitms = readtable("../PITMS_outputs/Yeast_avg_pitms_ColabFold_AF2Mv2_unpaired+paired.txt", opts);
 
load('../network_similarities/CH_L2_L3/results/CHA_L3_scores_net_perturbed_10percent_GSP_GSN_Yeast_DIP.mat');
load('../matrix/list_pairs_10percent_GSP_GSN_Yeast.mat');
load('../matrix/Yeast_DIP_net.mat');

%% Convert to output type
ints = table2cell(ints);
pitms = table2cell(pitms);

num_sims = 10;
INTS_eval_measures_mat = cell(num_sims, 1);
res_INTS_AUCmROC = zeros(num_sims, 1);
res_INTS_AUCROC = zeros(num_sims, 1);
res_INTS_AUCPR = zeros(num_sims, 1);
pairs_not_found = {};
count = 0;
for i = 1:num_sims
    pair_pos = L_pairs{i,1};
    pair_neg = L_pairs{i,2};
    assert(length(pair_pos) == length(pair_neg), 'Positive and negative sets have different size.');
    ints_pos = [];
    ints_neg = [];
    for j = 1:length(pair_pos)   
        idx = find(contains(ints(:,1), [map_lcc{pair_pos(j,1)}, '_', map_lcc{pair_pos(j,2)}]));
        if isempty(idx)
            idx = find(contains(ints(:,1), [map_lcc{pair_pos(j,2)}, '_', map_lcc{pair_pos(j,1)}]));
            if isempty(idx)
                count = count + 1;
                pairs_not_found = [pairs_not_found ; map_lcc{pair_pos(j,1)}, '_', map_lcc{pair_pos(j,2)}];
                disp([map_lcc{pair_pos(j,1)}, '_', map_lcc{pair_pos(j,2)}, ' not found (pos).']);
                continue;
            end
        end
        ints_pos = [ints_pos ; ints{idx,2}];
        
    end
    
    for j = 1:length(pair_neg) 
        idx = find(contains(ints(:,3), [map_lcc{pair_neg(j,1)}, '_', map_lcc{pair_neg(j,2)}]));
        if isempty(idx)
            idx = find(contains(ints(:,3), [map_lcc{pair_neg(j,2)}, '_', map_lcc{pair_neg(j,1)}]));
            if isempty(idx)
                count = count + 1;
                pairs_not_found = [pairs_not_found ; map_lcc{pair_neg(j,1)}, '_', map_lcc{pair_neg(j,2)}];
                disp([map_lcc{pair_neg(j,1)}, '_', map_lcc{pair_neg(j,2)}, ' not found (neg).']);
                continue;
            end
        end
        if isnan(ints{idx,4})
            disp([num2str(j), ' this is j']) 
            disp([map_lcc{pair_neg(j,1)}, '_', map_lcc{pair_neg(j,2)}, 'NAN?']);
        end
        ints_neg = [ints_neg ; ints{idx,4}];
        
    end
    labels = repelem([1, 0], [length(ints_pos) length(ints_neg)])';
    predicted_scores_INTS = [ints_pos; ints_neg];
    [measures_INTS, ~, ~] = prediction_evaluation(predicted_scores_INTS, labels);
    
    INTS_eval_measures_mat{i} = measures_INTS;
    res_INTS_AUCmROC(i) = measures_INTS.auc_mroc;
    res_INTS_AUCROC(i) = measures_INTS.auc_roc;
    res_INTS_AUCPR(i) = measures_INTS.auc_pr;
    
end

PITMS_eval_measures_mat = cell(num_sims, 1);
res_PITMS_AUCmROC = zeros(num_sims, 1);
res_PITMS_AUCROC = zeros(num_sims, 1);
res_PITMS_AUCPR = zeros(num_sims, 1);
count = 0;
for i = 1:num_sims
    pair_pos = L_pairs{i,1};
    pair_neg = L_pairs{i,2};
    assert(length(pair_pos) == length(pair_neg), 'Positive and negative sets have different size.');
    pitms_pos = [];
    pitms_neg = [];
    for j = 1:length(pair_pos)   
        idx = find(contains(pitms(:,1), [map_lcc{pair_pos(j,1)}, '_', map_lcc{pair_pos(j,2)}]));
        if isempty(idx)
            idx = find(contains(pitms(:,1), [map_lcc{pair_pos(j,2)}, '_', map_lcc{pair_pos(j,1)}]));
            if isempty(idx)
                count = count + 1;
                disp([map_lcc{pair_pos(j,1)}, '_', map_lcc{pair_pos(j,2)}, ' not found.']);
                continue;
            end
        end
        pitms_pos = [pitms_pos ; pitms{idx,2}];
        
    end
    
    for j = 1:length(pair_neg) 
        idx = find(contains(pitms(:,3), [map_lcc{pair_neg(j,1)}, '_', map_lcc{pair_neg(j,2)}]));
        if isempty(idx)
            idx = find(contains(pitms(:,3), [map_lcc{pair_neg(j,2)}, '_', map_lcc{pair_neg(j,1)}]));
            if isempty(idx)
                count = count + 1;
                disp([map_lcc{pair_neg(j,1)}, '_', map_lcc{pair_neg(j,2)}, ' not found.']);
                continue;
            end
        end
        pitms_neg = [pitms_neg ; pitms{idx,4}];
        
    end
    labels = repelem([1, 0], [length(pitms_pos) length(pitms_neg)])';
    predicted_scores_PITMS = [pitms_pos; pitms_neg];
    [measures_PITMS, ~, ~] = prediction_evaluation(predicted_scores_PITMS, labels);
    
    PITMS_eval_measures_mat{i} = measures_PITMS;
    res_PITMS_AUCmROC(i) = measures_PITMS.auc_mroc;
    res_PITMS_AUCROC(i) = measures_PITMS.auc_roc;
    res_PITMS_AUCPR(i) = measures_PITMS.auc_pr;
    
end

CH_methods = {'CH2-L2', 'CH3-L2', 'CH2-L3', 'CH3-L3'};
CH_eval_measures_mat = cell(num_sims, length(CH_methods));
res_CH_AUCmROC = zeros(num_sims, length(CH_methods));
res_CH_AUCROC = zeros(num_sims, length(CH_methods));
res_CH_AUCPR = zeros(num_sims, length(CH_methods));
for i = 1:num_sims
    
    S_tmp = table2array(removevars(S{i}, {'RA_L2_SPcorr', 'RA_L3_SPcorr', ...
                                          'CH1_L2_SPcorr', 'CH1_L3_SPcorr', ...
                                          'iLCL_L2_SPcorr', 'iLCL_L3_SPcorr'}));
    CH_mat = zeros(size(x_lcc, 1), size(x_lcc,2));
    ind = sub2ind([size(CH_mat, 1) size(CH_mat, 2)], S_tmp(:,1), S_tmp(:,2));
    P = L_pairs{i,1};
    N = L_pairs{i,2};
    labels = repelem([1, 0], size(P,1))';
    for j = 1:length(CH_methods)
        CH_mat(ind) = S_tmp(:,j+3);
        CH_mat_sym = CH_mat + CH_mat.';
        CH_pos_scores = cell(length(L_pairs{i}), 2);
        CH_neg_scores = cell(length(L_pairs{i}), 2);
        for k = 1:length(L_pairs{i})
            CH_pos_scores(k, :) = {[num2str(k) '_' num2str(P(k,1)) '_' num2str(P(k,2))] CH_mat_sym(P(k,1), P(k,2))};
            CH_neg_scores(k, :) = {[num2str(k) '_' num2str(N(k,1)) '_' num2str(N(k,2))] CH_mat_sym(N(k,1), N(k,2))};
            
        end
        predicted_scores_CH = [CH_pos_scores(:,2); CH_neg_scores(:,2)];
        [measures_CH, ~, ~] = prediction_evaluation(cell2mat(predicted_scores_CH), labels);
        CH_eval_measures_mat{i, j} = measures_CH;
        res_CH_AUCmROC(i, j) = measures_CH.auc_mroc;
        res_CH_AUCROC(i, j) = measures_CH.auc_roc;
        res_CH_AUCPR(i, j) = measures_CH.auc_pr;
    end
    
end

len = 4;
red = [1, 0, 0];
pink = [255, 192, 203]/255;
colors_p = [linspace(red(1),pink(1),len)', linspace(red(2),pink(2),len)', linspace(red(3),pink(3),len)'];
[y,i] = sort(mean(res_CH_AUCROC,1), 'descend');
methods = {CH_methods{i}, 'AFM-IS', 'AFM-piTMS'};
colors = {colors_p(1,:), colors_p(2,:), colors_p(3,:), colors_p(4,:), [0 0 0], [169 169 169]/255};
Y =  [mean(res_CH_AUCROC(:,i),1)'; mean(res_INTS_AUCROC); mean(res_PITMS_AUCROC)];
stderror_series = [std(res_CH_AUCROC(:,i),1)'/sqrt(size(res_CH_AUCROC, 1)); ...
    std(res_INTS_AUCROC)/sqrt(size(res_INTS_AUCROC, 1)); ...
    std(res_PITMS_AUCROC)/sqrt(size(res_PITMS_AUCROC, 1))];
[y,i] = sort(Y, 'descend');
subplot(2,1,1);
b = bar(y, 'FaceColor','flat');
b.CData(1,:) = colors{i(1)};
b.CData(2,:) = colors{i(2)};
b.CData(3,:) = colors{i(3)};
b.CData(4,:) = colors{i(4)};
b.CData(5,:) = colors{i(5)};
b.CData(6,:) = colors{i(6)};

ctr(1,:) = bsxfun(@plus, b(1).XData, b(1).XOffset');   
ydt(1,:) = b(1).YData;                                   
text(ctr(1,:), ydt(1,:) - 0.06, sprintfc('%.3f', ydt(1,:)), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',22, 'Color','k')
text(ctr(5), ydt(5) - 0.06, sprintfc('%.3f', ydt(5)), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',22, 'Color','w')

set(gca, 'FontSize', 13);
set(gca,'XTickLabel', methods(i), 'FontSize', 13);
ylim([0.5 1])
yline(max(y), '--', 'LineWidth', 3)
ylabel('AUC-ROC')

Pos = get(subplot(2,1,1),'Position');
P1=[ctr(5), ydt(5)];   
P2=[ctr(5), ydt(1)];         
Xlim=[0 6.5];
Ylim=[0.5 1];
hold on
X_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P1(1)-Xlim(1));
X_conv(2)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P2(1)-Xlim(1));
Y_conv(1)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P1(2)-Ylim(1));
Y_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P2(2)-Ylim(1));
xlim(Xlim)
ylim(Ylim)
annotation('arrow', X_conv, Y_conv, 'Color', [169 169 169]/255, 'LineWidth', 2.5)

perc = (ydt(1) - ydt(5)) / ydt(5) * 100;
text(ctr(5),ydt(1), sprintfc('+%.1f%%', perc), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',18, 'Color','r')

Pos = get(subplot(2,1,1),'Position');
P1=[ctr(6), ydt(6)];    
P2=[ctr(6), ydt(1)];     
Xlim=[0 6.5];
Ylim=[0.5 1];
hold on
X_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P1(1)-Xlim(1));
X_conv(2)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P2(1)-Xlim(1));
Y_conv(1)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P1(2)-Ylim(1));
Y_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P2(2)-Ylim(1));
xlim(Xlim)
ylim(Ylim)
annotation('arrow', X_conv, Y_conv, 'Color', [169 169 169]/255, 'LineWidth', 2.5)

perc = (ydt(1) - ydt(6)) / ydt(6) * 100;
text(ctr(6),ydt(1), sprintfc('+%.1f%%', perc), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',18, 'Color','r')

hold on
err = errorbar(1:6, Y(i), stderror_series(i), 'k', 'linestyle', 'none');
err.LineWidth = 2.5;
box off
hold off


[y,i] = sort(mean(res_CH_AUCPR,1), 'descend');
methods = {CH_methods{i}, 'AFM-IS', 'AFM-piTMS'};
Y =  [mean(res_CH_AUCPR(:,i),1)'; mean(res_INTS_AUCPR); mean(res_PITMS_AUCPR)];
stderror_series = [std(res_CH_AUCPR(:,i),1)'/sqrt(size(res_CH_AUCPR, 1)); ...
    std(res_INTS_AUCPR)/sqrt(size(res_INTS_AUCPR, 1)); ...
    std(res_PITMS_AUCPR)/sqrt(size(res_PITMS_AUCPR, 1))];
[y,i] = sort(Y, 'descend');
subplot(2,1,2);
b = bar(y, 'FaceColor','flat');
b.CData(1,:) = colors{i(1)};
b.CData(2,:) = colors{i(2)};
b.CData(3,:) = colors{i(3)};
b.CData(4,:) = colors{i(4)};
b.CData(5,:) = colors{i(5)};
b.CData(6,:) = colors{i(6)};

ctr(1,:) = bsxfun(@plus, b(1).XData, b(1).XOffset');    
ydt(1,:) = b(1).YData;                                   
text(ctr(1,:), ydt(1,:) - 0.06, sprintfc('%.3f', ydt(1,:)), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',22, 'Color','k')
text(ctr(5), ydt(5) - 0.06, sprintfc('%.3f', ydt(5)), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',22, 'Color','w')

set(gca, 'FontSize', 13);
set(gca,'XTickLabel', methods(i), 'FontSize', 13);
ylim([0.5 1])
yline(max(y), '--', 'LineWidth', 3)
xlabel('Methods') 
ylabel('AUC-PR')

Pos = get(subplot(2,1,2),'Position');
P1=[ctr(5), ydt(5)];  
P2=[ctr(5), ydt(1)];       
Xlim=[0 6.5];
Ylim=[0.5 1];
hold on
X_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P1(1)-Xlim(1));
X_conv(2)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P2(1)-Xlim(1));
Y_conv(1)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P1(2)-Ylim(1));
Y_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P2(2)-Ylim(1));
xlim(Xlim)
ylim(Ylim)
annotation('arrow', X_conv, Y_conv, 'Color', [169 169 169]/255, 'LineWidth', 2.5)

perc = (ydt(1) - ydt(5)) / ydt(5) * 100;
text(ctr(5),ydt(1), sprintfc('+%.1f%%', perc), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',18, 'Color','r')

Pos = get(subplot(2,1,2),'Position');
P1=[ctr(6), ydt(6)];     
P2=[ctr(6), ydt(1)];         
Xlim=[0 6.5];
Ylim=[0.5 1];
hold on
X_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P1(1)-Xlim(1));
X_conv(2)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P2(1)-Xlim(1));
Y_conv(1)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P1(2)-Ylim(1));
Y_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P2(2)-Ylim(1));
xlim(Xlim)
ylim(Ylim)
annotation('arrow', X_conv, Y_conv, 'Color', [169 169 169]/255, 'LineWidth', 2.5)

perc = (ydt(1) - ydt(6)) / ydt(6) * 100;
text(ctr(6),ydt(1), sprintfc('+%.1f%%', perc), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize',18, 'Color','r')

hold on
err = errorbar(1:6, Y(i), stderror_series(i), 'k', 'linestyle', 'none');
err.LineWidth = 2.5;
box off
hold off






