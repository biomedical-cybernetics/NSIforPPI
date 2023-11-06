function [matrices, percs_rem] = random_link_removal(x, percs, reps, method, filename)

%%% INPUT %%%
% x        - adjacency matrix
% percs    - vector of percentages of links to delete (example: [0.1, 0.2, 0.3])
% reps     - number of times that the removal process has to be repeated
% method   - removal method, the possible choices are described below
% filename - [optional] filename to store the output variables

%%% OUTPUT %%%
% matrices   - adjacency matrices obtained by the removal process
%              it is a cell of size [reps,length(percs)]
%              if it is a single matrix, it returns just the matrix (not a cell)
% percs_rem  - vector of percentages of links that have been removed
%              it is meaningful for methods that aim to keep the network connected
%              since the values could be lower than the desired input percentages

%%% METHODS %%%
% rand       - random deletion of a percentage of links
% rand_conn  - random deletion of a percentage of links keeping the network connected;
%              links to the leaves are iteratively excluded
%              since they would certainly disconnect the network
% DFS        - random deletion of a percentage of links
%              excluding a spanning tree to keep the network connected;
%              the spanning tree is obtained by the DFS algorithm
%              starting from a random root node;
%              starting from the same root node the algorithm returns always
%              the same spanning tree, it means that the spanning tree to preserve
%              is randomly chosen only among at most N spanning trees, not among all.
%              This gives a preference to certain links to be not removed,
%              therefore this algorithm has the disadvantage that the choice is biased.
%              In order to reduce the bias, for every repetition the node indexes
%              are randomly permutated before running the DFS and then restored.
%              This modification does not limit the retrievable spanning trees to N.
%              Since the DFS function runs a mex file, it is faster than the UST method.
% UST        - random deletion of a percentage of links
%              excluding a spanning tree to keep the network connected;
%              the uniform spanning tree is obtained by the Wilson'algorithm.
%              Since it is chosen randomly among all possible spanning trees 
%              with equal probability, the choice is not biased.
%              However, the method is slower than the DFS.

% mapping between input method names and function handles
methods = containers.Map();
methods('rand')      = @rand_link_rem_rand;
methods('rand_conn') = @rand_link_rem_rand_conn;
methods('DFS')       = @rand_link_rem_DFS;
methods('UST')       = @rand_link_rem_UST;

% check input arguments
narginchk(4, 5);
if any(percs < 0 | percs > 1)
    error('"percs" values must be in [0,1]')
end
if not(reps > 0 && mod(reps,1)==0)
    error('"reps" must be a positive integer')
end
if ~any(strcmp(method, keys(methods)))
    error(['the method "' method '" does not exist'])
end

% check input network
if ~issymmetric(double(x))
    x = max(x,x');
    warning('Input network has been transformed into undirected')
end
if any(x(:) ~= 0 & x(:) ~= 1)
    x(x ~= 0 & x ~= 1) = 1;
    warning('Input network has been transformed into unweighted')
end
if any(x(eye(size(x))==1))
    x(eye(size(x))==1) = 0;
    warning('Self-loops have been removed')
end
x = sparse(logical(x));

% if the method requires to preserve connectivity
% check that the input network is not already disconnected
if ~strcmp(method, 'rand')
    [ncc,~] = graphconncomp(x,'Directed',false);
    if ncc > 1
        error('Input network already disconnected, connectivity cannot be preserved')
    end
end

% initialization
matrices = cell(reps,length(percs));
percs_rem = zeros(reps,length(percs));

% convert function name to function handle
f = methods(method);

parfor_percent(reps);
parfor i = 1:reps
    
    [matrices(i,:), percs_rem(i,:)] = random_link_removal_rep_i(x, percs, f);
    perc = parfor_percent;
    if any(round(perc*reps) == round(reps/10.*(1:10)))
        display(['Percentage: ' num2str(round(perc*100))]);
    end
end
parfor_percent(0);

% the removed percentage is the same for all the repetitions
percs_rem = percs_rem(1,:);

% check removed percentages
for j = 1:length(percs)
    if percs(j) ~= percs_rem(j)
        warning(['Percentage to remove: ' num2str(percs(j)) '. Removable percentage: ' num2str(percs_rem(j)) '.'])
    end
end

% if there is a single matrix in output, do not return a cell
if reps==1 && length(percs)==1
    matrices = matrices{1,1};
end

% save output, if required
if nargin == 5
    save(filename, 'matrices', 'percs_rem', '-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Support Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function [matrices_i, percs_rem_i] = random_link_removal_rep_i(x, percs, f)

matrices_i = cell(1,length(percs));
percs_rem_i = zeros(1,length(percs));
for j = 1:length(percs)
    [x_rem, percs_rem_i(1,j)] = f(x, percs(j));
    matrices_i{1,j} = double(x_rem);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% method: rand

function [x, perc] = rand_link_rem_rand(x, perc)

% take upper triangle
x = triu(x,1);

[r,c] = find(x);            % find the existing links
w = [r c]; clear r c        % indexes of existing links
d = size(w,1);              % number of existing links

num = round(d*perc);        % number of links to delete

% randomly sample links to delete
idx = randsample(d, num);

% delete links
x(sub2ind(size(x), w(idx,1), w(idx,2))) = 0;

% return the matrix symmetric
x = max(x,x');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% method: rand_conn

function [x, perc] = rand_link_rem_rand_conn(x, perc)

% take upper triangle
x = triu(x,1);

d = sum(x(:));              % number of existing links
num = round(d*perc);        % number of links to delete

% check if the remaining links would be enough to keep the network connected
% for N nodes, at least N-1 links are required
if (d - num) < (size(x,1)-1)
    num = d - (size(x,1)-1);
    perc = num/d;
end

% iterative exclusion of leaves
xr = x;                                  % matrix indicating the removable links
leaves = (sum(xr,1)+sum(xr,2)'==1);      % find new leaves
while any(leaves)
    xr(leaves,:) = 0;                    % exclude leaves
    xr(:,leaves) = 0;
    leaves = (sum(xr,1)+sum(xr,2)'==1);  % find new leaves
end

[r,c] = find(xr);       % find the removable links
wr = [r c]; clear r c   % indexes of removable links
dr = size(wr,1);        % number of removable links

ncc = Inf;
while ncc > 1
    temp = x;
    
    % randomly sample links to delete
    idx = randsample(dr, num);
    
    % delete links
    temp(sub2ind(size(temp), wr(idx,1), wr(idx,2))) = 0;
    
    % compute the number of connected components
    % for undirected graphs the algorithm considers only the lower triangle
    [ncc,~] = graphconncomp(temp','Directed',false);
end

% return the matrix symmetric
x = temp; clear temp
x = max(x,x');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% method: DFS

function [x, perc] = rand_link_rem_DFS(x, perc)

d = sum(x(:))/2;                % number of existing links
num = round(d*perc);            % number of links to delete

% random permutation of the nodes
perm = randperm(size(x,1));
x_perm = x(perm,:);             % permute rows
x_perm = x_perm(:,perm);        % permute columns
map = [(1:size(x,1))' perm'];
map_rev = sortrows(map,2);
perm_rev = map_rev(:,1);

% choose randomly a root node
root = randsample(size(x,1), 1);

% compute a spanning tree using DFS
% the algorithm returns the predecessor for each node
% 0 is the predecessor of the root node
[~, pred, ~] = graphtraverse(x_perm, root, 'Directed', false, 'Method', 'DFS'); clear x_perm
tree = [(1:length(pred))' pred'];	% pairs node-predecessor
tree(root,:) = [];                  % remove fake pair for the root node
tree = sparse(tree(:,1), tree(:,2), ones(size(tree,1),1), size(x,1), size(x,2));    % create sparse directed tree
tree = max(tree,tree');				% make tree undirected

% restore the previous node indexes
tree = tree(perm_rev,:);
tree = tree(:,perm_rev);

% exclude the spanning tree
xr = x - tree;              % matrix indicating the removable links

[r,c] = find(triu(xr,1));   % find the removable links
wr = [r c]; clear r c       % indexes of removable links
dr = size(wr,1);            % number of removable links

if num > dr
    % if the desired percentage is not reachable
    % remove all the removable links
    x = x - xr;
    perc = dr/d;
else
    % otherwise randomly choose among the removable links
    idx = randsample(dr, num);
    x(sub2ind(size(x), wr(idx,1), wr(idx,2))) = 0;
    x(sub2ind(size(x), wr(idx,2), wr(idx,1))) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% method: UST

function [x, perc] = rand_link_rem_UST(x, perc)

d = sum(x(:))/2;            % number of existing links
num = round(d*perc);        % number of links to delete

% compute a uniform spanning tree
tree = unispantree_wilson(x);

% exclude the spanning tree
xr = x - tree;              % matrix indicating the removable links

[r,c] = find(triu(xr,1));   % find the removable links
wr = [r c]; clear r c       % indexes of removable links
dr = size(wr,1);            % number of removable links

if num > dr
    % if the desired percentage is not reachable
    % remove all the removable links
    x = x - xr;
    perc = dr/d;
else
    % otherwise randomly choose among the removable links
    idx = randsample(dr, num);
    x(sub2ind(size(x), wr(idx,1), wr(idx,2))) = 0;
    x(sub2ind(size(x), wr(idx,2), wr(idx,1))) = 0; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tree = unispantree_wilson(x)

% The function computes a uniform spanning tree, a spanning tree
% chosen randomly among all possible spanning trees with equal probability.
% It is based on the algorithm RandomTreeWithRoot proposed in:
% Wilson D.B. (1996), "Generating Random Spanning Trees More Quickly than the Cover Time"
% It maintains a current tree, which initially consists of just the root.
% While there are remaining vertices not in the tree, the algorithm
% does a random walk from one such vertex, erasing cycles as they are created,
% until the walk hits the current tree. Then the cycle erased trajectory
% gets added to the current tree.
% The computational complexity is O(tau), where tau is the mean hitting time,
% the expected time it takes for a random walk between two random vertices.
%
% In this implementation:
% - the input network is assumed to be unweighted, undirected and connected
% - the root is randomly chosen among the nodes with probability
%   proportional to their degree (following the suggestion in the paper to
%   choose as root a random endpoint of a random edge to make the algorithm
%   more efficient in case of undirected graphs)

%%% INPUT %%%
% x - adjacency matrix of the network

%%% OUTPUT %%%
% tree - adjacency matrix of the uniform spanning tree

% initialization
n = size(x,1);              % number of nodes
in_tree = false(n,1);       % nodes currently in the tree
next = zeros(n,1);          % successors of nodes
deg = full(sum(x,1));       % node degree

% compute adjacency list and node degrees
adj_list = cell(n,1);
for i = 1:n
    adj_list{i} = find(x(i,:));
end

% randomly choose a root with probability proportional to the node degree
r = randsample(n,1,true,deg);

% add the root to the tree
in_tree(r) = true;

% random node sequence
rand_seq = randperm(n);

for i = 1:n
    % random walk from a random node until hitting the current tree
    % erasing cycles as they are created
    v = rand_seq(i);
    while ~in_tree(v)
        next(v) = adj_list{v}(randsample(deg(v),1));
        v = next(v);
    end
    % add the cycle-erased random walk to the current tree
    v = rand_seq(i);
    while ~in_tree(v)
        in_tree(v) = true;
        v = next(v);
    end
end

% build the tree adjacency matrix
tree = [(1:n)' next];               % pairs node-successor
tree(r,:) = [];                     % remove fake pair for the root node
tree = sparse(tree(:,1), tree(:,2), ones(size(tree,1),1), n, n);    % create sparse unweighted and directed tree
tree = max(tree,tree');				% make tree undirected

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function percent = parfor_percent(N)
%
%   Function for monitoring the percentage of the parfor progress.
%   It works by creating a file called parfor_progress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   1) Call PARFOR_PERCENT(N) before the parfor loop, N is the number of iterations.
%   It initializes the text file.
%
%   2) Call PARFOR_PERCENT inside the parfor loop to obtain the current percentage.
%   It updates and read the text file.
%
%   3) Call PARFOR_PERCENT(0) after the parfor loop.
%   It deletes the text file.
%
%   Example:
%
%      N = 100;
%      parfor_percent(N);
%      parfor i = 1:N
%         pause(rand); % Replace with real code
%         percent = parfor_percent;
%         percent
%      end
%      parfor_percent(0);
%
%   The original function is provided by Jeremy Scheff
%   In this modified version the printing of the progress bar is avoided

narginchk(0, 1);

if nargin == 0
    N = -1;
end

percent = 0;

if N > 0
    % Initialization
    f = fopen('parfor_progress.txt', 'w');
    if f < 0
        error('Do you have write permissions for %s?', pwd);
    end
    fprintf(f, '%d\n', N);
    fclose(f);
elseif N == 0
    % Finalization
    delete('parfor_progress.txt');
    percent = 100;
else
    % Updating
    if ~exist('parfor_progress.txt', 'file')
        error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
    end
    
    f = fopen('parfor_progress.txt', 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen('parfor_progress.txt', 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1);
end