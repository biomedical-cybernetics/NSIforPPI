function [x_rem, perc_rem] = random_link_removal(x, perc, method, rand_conn_maxiters, check_conn, show_warnings)

%%% INPUT %%%
% x      - adjacency matrix
% perc   - percentage of links to delete
% method - removal method, the possible choices are described below
% rand_conn_maxiters - maximum number of attempts for the 'rand_conn' method
%                      [optional, if empty or not given the default is Inf]
% check_conn - 1 or 0 to check if the network is connected or not
%              [optional, if empty or not given the default is 1]
% show_warnings - 1 or 0 to show the warnings or not
%              [optional, if empty or not given the default is 1]

%%% OUTPUT %%%
% x_rem    - adjacency matrix obtained by the removal process
% perc_rem - percentage of links that have been removed
%            it is meaningful for methods that aim to keep the network connected
%            since the value could be lower than the desired input percentage

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
% mindeg1    - random deletion of a percentage of links ensuring minimum degree 1;
%              for each node, one link is preserved at random.

% check input arguments
narginchk(3,6);
validateattributes(x, {'numeric'}, {'square','binary'});
if ~issymmetric(x)
    error('Input network must be undirected')
end
if any(x(speye(size(x))==1))
    error('Input network must be zero-diagonal')
end
x = sparse(x);
validateattributes(perc, {'numeric'}, {'scalar','>=',0,'<=',1});
validateattributes(method, {'char'}, {});
if ~any(strcmp(method, {'rand','rand_conn','DFS','UST','mindeg1'}))
    error(['the method "' method '" does not exist'])
end
if ~exist('rand_conn_maxiters', 'var') || isempty(rand_conn_maxiters)
    rand_conn_maxiters = Inf;
else
    validateattributes(rand_conn_maxiters, {'numeric'}, {'scalar','integer','positive'});
end
if ~exist('check_conn', 'var') || isempty(check_conn)
    check_conn = 1;
else
    validateattributes(check_conn, {'numeric'}, {'scalar','binary'});
end
if ~exist('show_warnings', 'var') || isempty(show_warnings)
    show_warnings = 1;
else
    validateattributes(show_warnings, {'numeric'}, {'scalar','binary'});
end

% check that the input network is connected
if check_conn
    [ncc,~] = graphconncomp(x, 'Directed', false);
    if ncc > 1
        error('Input network already disconnected, connectivity cannot be preserved')
    end
end

% convert function name to function handle
if strcmp(method, 'rand')
    [x_rem, perc_rem] = rand_link_rem_rand(x, perc);
elseif strcmp(method, 'rand_conn')
    [x_rem, perc_rem] = rand_link_rem_rand_conn(x, perc, rand_conn_maxiters);
elseif strcmp(method, 'DFS')
    [x_rem, perc_rem] = rand_link_rem_DFS(x, perc);
elseif strcmp(method, 'UST')
    [x_rem, perc_rem] = rand_link_rem_UST(x, perc);
elseif strcmp(method, 'mindeg1')
    [x_rem, perc_rem] = rand_link_rem_mindeg1(x, perc);
end

% warnings
if show_warnings
    if perc_rem < perc
        warning(['Percentage to remove: ' num2str(perc) '. Removable percentage: ' num2str(perc_rem) '.'])
    end
    if isempty(x_rem)
        warning('rand_conn: maximum iterations reached.')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Support Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [x, perc] = rand_link_rem_rand_conn(x, perc, rand_conn_maxiters)

% take upper triangle
x = triu(x,1);

d = sum(sum(x));            % number of existing links
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
iter = 1;
while ncc > 1 && iter <= rand_conn_maxiters
    temp = x;
    
    % randomly sample links to delete
    idx = randsample(dr, num);
    
    % delete links
    temp(sub2ind(size(temp), wr(idx,1), wr(idx,2))) = 0;
    
    % compute the number of connected components
    % for undirected graphs the algorithm considers only the lower triangle
    [ncc,~] = graphconncomp(temp','Directed',false);
    
    iter = iter + 1;
end

if ncc == 1
    % return the matrix symmetric
    x = temp; clear temp
    x = max(x,x');
else
    x = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% method: DFS

function [x, perc] = rand_link_rem_DFS(x, perc)

d = sum(sum(x))/2;              % number of existing links
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

d = sum(sum(x))/2;          % number of existing links
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

% method: mindeg1

function [x, perc] = rand_link_rem_mindeg1(x, perc)

[~, idx] = sort(sum(x,2));
mask = zeros(length(x),1);
to_keep = [];
for i = 1:length(x)
    n1 = idx(i);
    if mask(n1)==1
        continue;
    end
    A = find(x(n1,:));
    n2 = A(randsample(length(A), 1));
    to_keep = [to_keep; [min(n1,n2),max(n1,n2)]];
    mask(n1) = 1;
    mask(n2) = 1;
end

% take upper triangle
x = triu(x,1);

[r,c] = find(x);            % find the existing links
w = [r c]; clear r c        % indexes of existing links
d = size(w,1);              % number of existing links
num = round(d*perc);        % number of links to delete

w = setdiff(w, to_keep, 'rows');
if num <= size(w,1)
    % randomly sample links to delete
    idx = randsample(size(w,1), num);
else
    idx = 1:size(w,1);
    perc = size(w,1)/d;
end

% delete links
x(sub2ind(size(x), w(idx,1), w(idx,2))) = 0;

% return the matrix symmetric
x = max(x,x');
