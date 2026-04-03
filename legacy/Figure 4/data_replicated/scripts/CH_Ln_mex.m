function [w, scores_Ln, scores_RA_Ln, scores_CH1_Ln, scores_CH2_Ln, scores_CH3_Ln, scores_CH4_Ln] = CH_Ln_mex(x, Ln, w, par)

%%% INPUT %%%
% x - adjacency matrix of the network
%     the network is considered unweighted, undirected and zero-diagonal
% Ln - maximum length of paths
% w - [optional] two-columns matrix (id1,id2) indicating the links for which the score should be calculated
%     if not given or empty, the scores for all the missing links are computed
% par - [optional] 1 or 0 to indicate whether the function should use parallel computation
%       if not given, parallel computation is used
%
%%% OUTPUT %%%
% w - two-columns matrix (id1,id2) indicating the links for which the scores have been calculated
% scores_Ln     - (Ln-1)-columns matrix containing the scores for L2,L3,...,Ln
% scores_RA_Ln  - (Ln-1)-columns matrix containing the scores for RA_L2,RA_L3,...,RA_Ln
% scores_CH1_Ln - (Ln-1)-columns matrix containing the scores for CH1_L2,CH1_L3,...,CH1_Ln
% scores_CH2_Ln - (Ln-1)-columns matrix containing the scores for CH2_L2,CH2_L3,...,CH2_Ln
% scores_CH3_Ln - (Ln-1)-columns matrix containing the scores for CH3_L2,CH3_L3,...,CH3_Ln
% scores_CH4_Ln - (Ln-1)-columns matrix containing the scores for CH4_L2,CH4_L3,...,CH4_Ln

% check input
narginchk(2,4)
validateattributes(x, {'numeric'}, {'square','binary'});
x = sparse(max(x,x'));
x(speye(size(x))==1) = 0;
validateattributes(Ln, {'numeric'}, {'scalar','integer','>=',2});
if ~exist('w', 'var') || isempty(w)
    [e1,e2] = find(triu(x==0,1));
    w = [e1,e2];
else
    validateattributes(w, {'numeric'}, {'ncols',2,'integer','positive'});
    e1 = w(:,1); e2 = w(:,2);
end
if ~exist('par', 'var') || isempty(par)
    par = Inf;
else
    validateattributes(par, {'numeric'}, {'scalar','integer','>=',0,'<=',1});
    if par == 1
        par = Inf;
    end
end

% initialization
m = size(w,1);
scores_Ln = zeros(m,Ln-1);
scores_RA_Ln = zeros(m,Ln-1);
scores_CH1_Ln = zeros(m,Ln-1);
scores_CH2_Ln = zeros(m,Ln-1);
scores_CH3_Ln = zeros(m,Ln-1);
scores_CH4_Ln = zeros(m,Ln-1);
deg = full(sum(x,1));

% compute scores
parfor (i = 1:m, par)
    [scores_Ln(i,:), scores_RA_Ln(i,:), scores_CH1_Ln(i,:), scores_CH2_Ln(i,:), scores_CH3_Ln(i,:), scores_CH4_Ln(i,:)] = compute_scores_i(x, e1(i), e2(i), Ln, deg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [scores_Ln, scores_RA_Ln, scores_CH1_Ln, scores_CH2_Ln, scores_CH3_Ln, scores_CH4_Ln] = compute_scores_i(x, s, t, Ln, deg)

paths = find_paths_mex(x, s, t, Ln);

scores_Ln = zeros(1,Ln-1);
scores_RA_Ln = zeros(1,Ln-1);
scores_CH1_Ln = zeros(1,Ln-1);
scores_CH2_Ln = zeros(1,Ln-1);
scores_CH3_Ln = zeros(1,Ln-1);
scores_CH4_Ln = zeros(1,Ln-1);
for l = 1:Ln-1
    if ~isempty(paths{l})
        paths_l = paths{l};
        paths_l_size = size(paths_l);

        if l == 1
            paths_l = reshape(paths_l, 1, numel(paths_l));
            degi_l = full(sum(x(paths_l,paths_l),1));
        else
            idx = zeros(2,paths_l_size(2)*(l-1));
            for j = 1:l-1
               idx(:,(j-1)*paths_l_size(2)+1:j*paths_l_size(2)) = paths_l(j:j+1,:);
            end
            xi = sparse([idx(1,:) idx(2,:)],[idx(2,:) idx(1,:)],1,length(x),length(x));
            degi_l = full(sum(xi>0,1));
            paths_l = reshape(paths_l, 1, numel(paths_l));
            degi_l = degi_l(paths_l);
        end
        deg_l = deg(paths_l);
        dege_l = deg_l - degi_l - full(x(s,paths_l)) - full(x(t,paths_l));
        
        deg_l = reshape(deg_l, paths_l_size);
        degi_l = reshape(degi_l, paths_l_size);
        dege_l = reshape(dege_l, paths_l_size);
        
        scores_Ln(l) = paths_l_size(2);
        scores_RA_Ln(l) = sum(1./geomean(deg_l,1),2);
        scores_CH1_Ln(l) = sum(geomean(degi_l,1)./geomean(deg_l,1),2);
        scores_CH2_Ln(l) = sum(geomean(1+degi_l,1)./geomean(1+dege_l,1),2);
        scores_CH3_Ln(l) = sum(1./geomean(1+dege_l,1),2);
        scores_CH4_Ln(l) = sum(geomean(degi_l,1)./geomean(1+degi_l,1),2);
    end
end
