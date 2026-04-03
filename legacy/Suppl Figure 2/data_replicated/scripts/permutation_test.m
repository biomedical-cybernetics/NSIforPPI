function p = permutation_test(x1, x2, func_name, tail_type, iters, rstr_flag)

%%% INPUT %%%
% x1, x2 - numerical vectors
% func_name - name of the function to evaluate the statistic for each vector (examples: 'mean', 'median')
% tail_type - choose between:
%             'both'  -> test for statistic in x1 greater or lower than in x2
%             'right' -> test for statistic in x1 greater than in x2
%             'left'  -> test for statistic in x1 lower than in x2
%             (if not given or empty, default is 'both')
% iters - number of iterations (if not given or empty, default = 1000)
% rstr_flag - 1 or 0 to indicate if a random stream generator should be used or not
%             for reproducibility of the results (if not given or empty, default = 0)

%%% OUTPUT %%%
% p - empirical p-value obtained using the formula discussed in:
%     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC379178/

% check input
narginchk(3,6)
validateattributes(x1, {'numeric'}, {'vector','finite'});
validateattributes(x2, {'numeric'}, {'vector','finite'});
validateattributes(func_name, {'char'}, {});
func = str2func(func_name);
if ~exist('tail_type', 'var') || isempty(tail_type)
    tail_type = 'both';
elseif ~any(strcmp(tail_type,{'right','left','both'}))
    error('Invalid value for ''tail_type''.')
end
if ~exist('iters', 'var') || isempty(iters)
    iters = 1000;
else
    validateattributes(iters, {'numeric'}, {'scalar','integer','positive'});
end
if ~exist('rstr_flag', 'var') || isempty(rstr_flag)
    rstr_flag = 0;
else
    validateattributes(rstr_flag, {'numeric'}, {'scalar','binary'});
    if rstr_flag
        rstr = RandStream('mt19937ar','Seed',1);
    end
end
if isrow(x1)
    x1 = x1';
end
if isrow(x2)
    x2 = x2';
end

% observed difference
diff_obs = func(x1) - func(x2);

% random distribution of differences
x = [x1; x2];
n = length(x);
n1 = length(x1);
diff_rand = zeros(1,iters);
for i = 1:iters
    if rstr_flag
        rp = randperm(rstr, n);
    else
        rp = randperm(n);
    end
    xr1 = x(rp(1:n1));
    xr2 = x(rp(n1+1:end));
    diff_rand(i) = func(xr1) - func(xr2);
end

% p-value
if strcmp(tail_type, 'right')
    p = (sum(diff_rand >= diff_obs) + 1) / (iters + 1);
elseif strcmp(tail_type, 'left')
    p = (sum(diff_rand <= diff_obs) + 1) / (iters + 1);
elseif strcmp(tail_type, 'both')
    p = (sum(abs(diff_rand) >= abs(diff_obs)) + 1) / (iters + 1);
end
