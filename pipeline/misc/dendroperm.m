function perm = dendroperm(Z,varargin)
%DENDROGRAM Generate dendrogram plot.
%   DENDROGRAM(Z) generates a dendrogram plot of the hierarchical binary
%   cluster tree represented by Z.  Z is an (M-1)-by-3 matrix, generated by
%   the LINKAGE function, where M is the number of objects in the original
%   dataset.
%
%   A dendrogram consists of many U-shaped lines connecting objects in a
%   hierarchical tree.  The height of each U represents the distance between
%   the two objects being connected.  If there were 30 or fewer data points in
%   the original dataset, each leaf in the dendrogram corresponds to one data
%   point.  If there were more than 30 data points, the complete tree can look
%   crowded, and DENDROGRAM collapses lower branches as necessary, so that
%   some leaves in the plot correspond to more than one data point.
%
%   DENDROGRAM(Z,P) generates a dendrogram with no more than P leaf nodes, by
%   collapsing lower branches of the tree.  To display the complete tree, set
%   P = 0.
%
%   H = DENDROGRAM(...) returns a vector of line handles.
%
%   [H,T] = DENDROGRAM(...) generates a dendrogram and returns T, a vector
%   of size M that contains the leaf node number for each object in the
%   original dataset.  T is useful when P is less than the total number of
%   objects, so some leaf nodes in the display correspond to multiple
%   objects.  For example, to find out which objects are contained in leaf
%   node k of the dendrogram, use find(T==k). When there are fewer than P
%   objects in the original data, all objects are displayed in the
%   dendrogram.  In this case, T is the identity map, i.e., T = (1:M)',
%   where each node contains only a single object.
%
%   [H,T,PERM] = DENDROGRAM(...) generates a dendrogram and returns the
%   permutation vector of the node labels of the leaves of the dendrogram.
%   PERM is ordered from left to right on a horizontal dendrogram and
%   bottom to top for a vertical dendrogram.
%
%   [...] = DENDROGRAM(...,'COLORTHRESHOLD',T) assigns a unique color to
%   each group of nodes within the dendrogram whose linkage is less than
%   the scalar value T where T is in the range 0 < T < max(Z(:,3)). If T is
%   less than or equal to zero or if T is greater than the maximum linkage
%   then the dendrogram will be drawn using only one color. T can also be
%   set to 'default' in which case the threshold will be set to 70% of the
%   maximum linkage i.e. 0.7 * max(Z(:,3)).
%
%   [...] = DENDROGRAM(...,'ORIENTATION',ORIENT) will orient the dendrogram
%   within the figure window. Options are:
%
%      'top'      --- top to bottom (default)
%      'bottom'   --- bottom to top
%      'left'     --- left to right
%      'right'    --- right to left
%
%   [...] = DENDROGRAM(...,'LABELS',S) accepts a character array or cell array
%   of strings S with one label for each observation.  Any leaves in the tree
%   containing a single observation are labeled with that observation's label.
%
%   Example:
%      X = rand(100,2);
%      Y = pdist(X,'cityblock');
%      Z = linkage(Y,'average');
%      [H, T] = dendrogram(Z);
%
%   See also LINKAGE, PDIST, CLUSTER, CLUSTERDATA, COPHENET, INCONSISTENT,
%   SILHOUETTE.

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.4 $

m = size(Z,1)+1;
if nargin < 2
    p = 30;
end

if nargin == 2
    p = varargin{1};
end

color = false;
orientation = 't';
obslabels = [];
threshold = 0.7 * max(Z(:,3));
leafOrder = [];
horz = false;

if nargin > 2
    if isnumeric(varargin{1})
        p = varargin{1};
        offset = 1;
    else
        p = 30;
        offset = 0;
    end
    
    pnames = {'orientation' 'colorthreshold' 'labels'  'reorder'};
    dflts =  {orientation   'default'        obslabels leafOrder};
    [orientation,threshold,obslabels,leafOrder,setFlag] =  ...
        internal.stats.parseArgs(pnames, dflts, varargin{1+offset:end});
    
    if ~isempty(orientation) && ischar(orientation)
        orientation = lower(orientation(1));
    else
        orientation = 0;    % bad value
    end

    if ~ismember(orientation,{'t','b','r','l'})
        orientation = 't';
        warning(message('stats:dendrogram:BadOrientation'));
    end
    horz = ismember(orientation,{'r','l'});
    
    if ~isempty(obslabels)
        if ischar(obslabels)
            obslabels = cellstr(obslabels);
        elseif ~iscellstr(obslabels)
            error(message('stats:dendrogram:BadLabels'));
        end
        if ~isvector(obslabels) || numel(obslabels)~=m
            error(message('stats:dendrogram:InputSizeMismatch'));
        end
        obslabels = obslabels(:);
    end
    if ~isempty(leafOrder) && (~isvector(leafOrder) || numel(leafOrder)~=m)
        error(message('stats:dendrogram:BadLeafOrder'));
    end
end
% For each node currently labeled m+k, replace its index by
% min(i,j) where i and j are the nodes under node m+k.
Z = transz(Z);
T = (1:m)';

% If there are more than p nodes, the dendrogram looks crowded.
% The following code will make the last p link nodes into leaf nodes,
% and only these p nodes will be visible.
if (m > p) && (p ~= 0)
    
    Y = Z((m-p+1):end,:);         % get the last nodes
    
    R = unique(Y(:,1:2));
    Rlp = R(R<=p);
    Rgp = R(R>p);
    W(Rlp) = Rlp;                 % use current node number if <=p
    W(Rgp) = setdiff(1:p, Rlp);   % otherwise get unused numbers <=p
    W = W(:);
    T(R) = W(R);
    
    % Assign each leaf in the original tree to one of the new node numbers
    for i = 1:p
        c = R(i);
        T = clusternum(Z,T,W(c),c,m-p+1,0); % assign to its leaves.
    end
    
    % Create new, smaller tree Z with new node numbering
    Y(:,1) = W(Y(:,1));
    Y(:,2) = W(Y(:,2));
    Z = Y;
    
    m = p; % reset the number of node to be 30 (row number = 29).
end

A = zeros(4,m-1);
B = A;
n = m;
X = 1:n;
Y = zeros(n,1);
r = Y;

% arrange Z into W so that there will be no crossing in the dendrogram.
W = zeros(size(Z));
W(1,:) = Z(1,:);

nsw = zeros(n,1); rsw = nsw;
nsw(Z(1,1:2)) = 1; rsw(1) = 1;
k = 2; s = 2;

while (k < n)
    i = s;
    while rsw(i) || ~any(nsw(Z(i,1:2)))
        if rsw(i) && i == s
            s = s+1;
        end
        i = i+1;
    end
    
    W(k,:) = Z(i,:);
    nsw(Z(i,1:2)) = 1;
    rsw(i) = 1;
    if s == i
        s = s+1;
    end
    k = k+1;
end

g = 1;
for k = 1:m-1 % initialize X
    i = W(k,1);
    if ~r(i),
        X(i) = g;
        g = g+1;
        r(i) = 1;
    end
    i = W(k,2);
    if ~r(i),
        X(i) = g;
        g = g+1;
        r(i) = 1;
    end
end

% if a leaf order is specified use
if ~isempty(leafOrder)
    [dummy, X] = sort(leafOrder); %#ok
end

[u,perm]=sort(X);   %#ok perm is the third output value


% ---------------------------------------
function T = clusternum(X, T, c, k, m, d)
% assign leaves under cluster c to c.

d = d+1;
n = m; flag = 0;
while n > 1
    n = n-1;
    if X(n,1) == k % node k is not a leave, it has subtrees
        T = clusternum(X, T, c, k, n,d); % trace back left subtree
        T = clusternum(X, T, c, X(n,2), n,d);
        flag = 1; break;
    end
end

if flag == 0 && d ~= 1 % row m is leaf node.
    T(X(m,1)) = c;
    T(X(m,2)) = c;
end
% ---------------------------------------
function T = colorcluster(X, T, k, m)
% find local clustering

n = m;
while n > 1
    n = n-1;
    if X(n,1) == k % node k is not a leave, it has subtrees
        T = colorcluster(X, T, k, n); % trace back left subtree
        T = colorcluster(X, T, X(n,2), n);
        break;
    end
end
T(m) = 1;
% ---------------------------------------
function Z = transz(Z)
%TRANSZ Translate output of LINKAGE into another format.
%   This is a helper function used by DENDROGRAM and COPHENET.

%   In LINKAGE, when a new cluster is formed from cluster i & j, it is
%   easier for the latter computation to name the newly formed cluster
%   min(i,j). However, this definition makes it hard to understand
%   the linkage information. We choose to give the newly formed
%   cluster a cluster index M+k, where M is the number of original
%   observation, and k means that this new cluster is the kth cluster
%   to be formmed. This helper function converts the M+k indexing into
%   min(i,j) indexing.

m = size(Z,1)+1;

for i = 1:(m-1)
    if Z(i,1) > m
        Z(i,1) = traceback(Z,Z(i,1));
    end
    if Z(i,2) > m
        Z(i,2) = traceback(Z,Z(i,2));
    end
    if Z(i,1) > Z(i,2)
        Z(i,1:2) = Z(i,[2 1]);
    end
end


function a = traceback(Z,b)

m = size(Z,1)+1;

if Z(b-m,1) > m
    a = traceback(Z,Z(b-m,1));
else
    a = Z(b-m,1);
end
if Z(b-m,2) > m
    c = traceback(Z,Z(b-m,2));
else
    c = Z(b-m,2);
end

a = min(a,c);
