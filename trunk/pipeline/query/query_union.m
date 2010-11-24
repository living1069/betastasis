
% QUERY_UNION     Merge samples from two or more query sets into one merged set
% 
%    MERGED = QUERY_UNION(A, B) merges the samples from query sets A and B and
%    returns the merged query set MERGED. Note that the function only works
%    with query sets, not with realized datasets. Any number of query sets
%    can be given as input to the function. For instance, QUERY_UNION(A, B, C)
%    would merge samples from three query sets.
%
%    Query sets can only be merged if they are of the same data type. For
%    instance, you can merge two gene expression query sets, but you cannot
%    merge a gene expression query set with a microRNA expression query set.
%    As an additional constraint, raw microarray probe intensity query sets can
%    only be merged if all samples come from the same platform.
%
%    If the input query sets contain samples with identical resource IDs,
%    duplicate samples are discarded so that only one instance of every
%    resource ID is present in the merged query set.

% Author: Matti Annala <matti.annala@tut.fi>

function qset = query_union(varargin)

if length(varargin) < 2
	error 'query_union() requires at least two input arguments.';
end

type = varargin{1}.Type;

for k = 1:length(varargin)
	if ~isfield(varargin{k}, 'Resource')
		error('Argument #%d of query_union() is not a query set.', k);
	end
	if ~strcmp(type, varargin{k}.Type)
		error('Arguments of query_union() must have the same query set type.');
	end
	if strcmp(type, 'Microarray probe intensities')
		if ~strcmp(varargin{1}.Platform{1}, varargin{k}.Platform{1})
			error(['Microarray probe intensity query sets can only be ' ...
			       'merged if they come from the same platform.']);
		end
	end
end

qset = cat_structs(varargin{:});

keep = true(length(qset.Resource), 1);
for k = 1:length(qset.Resource)
	if keep(k) == false, continue, end
	dups = find(strcmp(qset.Resource{k}, qset.Resource));
	dups = dups(2:end);
	keep(dups) = false;
end

if ~all(keep)
	fprintf(1, ['WARNING: Discarded %d duplicated samples ' ...
		'(identical resource IDs).\n'], sum(~keep));
end

qset = filter_struct(qset, keep);

