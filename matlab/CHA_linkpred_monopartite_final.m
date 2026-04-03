function scores = CHA_linkpred_monopartite_final(x, methods, cores)
narginchk(2,3)
validateattributes(x, {'numeric'}, {'square','binary'});
x = sparse(x);
if ~issymmetric(x)
    error('The input matrix must be symmetric.')
end
if any(x(speye(size(x))==1))
    error('The input matrix must be zero-diagonal.')
end
if ~exist('methods', 'var') || isempty(methods)
    methods = {'CH2_L2', 'CH3_L2', 'CH2_L3', 'CH3_L3'};
else
    validateattributes(methods, {'cell'}, {});
    if length(methods) > length(unique(methods))
        error('The variable ''methods'' should not contain duplicates.')
    end
end
if ~exist('cores', 'var') || isempty(cores)
    cores = Inf;
else
    validateattributes(cores, {'numeric'}, {'scalar','integer','positive'});
end

M = length(methods);
L = NaN(M,1);
models = cell(M,1);
for m = 1:M
    temp = strsplit(methods{m},'_L');
    L(m) = str2double(temp{2});
    models{m} = temp{1};
end
L = unique(L);
models_all = {'RA','CH1','CH2','CH3', 'iLCL', 'CH3.1'};
models = find(ismember(models_all,models))-1;
if isinf(cores)
    scores = CH_scores_mex_final(x, L, models);
else
    scores = CH_scores_mex_final(x, L, models, cores);
end
S = cell(length(methods),1);
for i1 = 1:length(L)
    for i2 = 1:length(models)
        method = sprintf('%s_L%d', models_all{models(i2)+1}, L(i1));
        m = find(strcmp(method,methods),1);
        if ~isempty(m)
            S{m} = scores{i1,i2};
            scores{i1,i2} = [];
        end
    end
end

if cores == 1 || M == 1
    cores = 0;
end
[e1,e2] = find(triu(true(size(x)),1));

% Save raw RA_L3 scores before SPcorr subranking
idx_ra_l3 = find(strcmp(methods, 'RA_L3'), 1);
if ~isempty(idx_ra_l3)
    RA_L3_Istvan_matrix = S{idx_ra_l3};
    scores_RA_L3_Istvan = RA_L3_Istvan_matrix(sub2ind(size(x),e1,e2));
end

scores = zeros(length(e1),M);
parfor (m = 1:M, cores)
    tmp = compute_Pcorr_and_rank_scores(x, S{m}, e1, e2);
    tmp = tiedrank(tmp);
    scores(:,m) = tmp;
end

colnames = {'node1', 'node2'};
for m = 1:M
    if m == idx_ra_l3
        colnames{end+1} = 'RA_L3';
    end
    colnames{end+1} = [methods{m} '_SPcorr'];
end

if ~isempty(idx_ra_l3)
    scores = array2table([e1 e2 scores(:, 1:idx_ra_l3-1) scores_RA_L3_Istvan scores(:, idx_ra_l3:end)],'VariableNames',colnames);
else
    scores = array2table([e1 e2 scores],'VariableNames',colnames);
end


function scores = compute_Pcorr_and_rank_scores(x, S, e1, e2)
s = x .* abs(S - min(S(S>0))-max(S(:)));
G = graph(sparse(s));
s = distances(G);
s = replace_inf_distances(s);
s = corrcoef(s);
scores = [S(sub2ind(size(x),e1,e2)) s(sub2ind(size(x),e1,e2))];
[~,~,scores] = unique(scores, 'sorted', 'rows');
