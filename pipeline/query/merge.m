
% MERGE     Merge samples from two or more query sets into one merged set
% 
%    MERGED = MERGE(A, B) merges the samples from query sets A and B and
%    returns the merged query set MERGED. Note that the function only works
%    with query sets, not with realized datasets. Any number of query sets
%    can be given as input to the function. For instance, MERGE(A, B, C)
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

function merged = merge(varargin)

merged = struct;

if length(varargin) < 2
	error 'merge() requires at least two input arguments.';
end

type = varargin{1}.meta.type;
for k = 2:length(varargin)
	if ~strcmp(type, varargin{k}.meta.type)
		error 'Datasets must be of the same type.';
	end
end

% Figure out what kind of a dataset we are dealing with and perform
% necessary tasks.
if rx(type, 'gene expression')
	for d = 1:length(varargin), labels{d} = varargin{d}.rows.gene_symbol; end
	for d = 1:length(varargin), matrices{d} = varargin{d}.mean; end
	[merged.rows.gene_symbol, merged.mean] = join_matrices(labels, matrices);
elseif rx(type, 'fusion genes')
	for d = 1:length(varargin), fusions{d} = varargin{d}.fusions; end
	for d = 1:length(varargin), ptl{d} = varargin{d}.paired_tag_length; end
	merged.fusions = cat(2, fusions{:});
	merged.paired_tag_length = cat(2, ptl{:});
else
	error 'Dataset is of unsupported type.';
end

metas = cell(1, length(varargin));
for d = 1:length(varargin), metas{d} = varargin{d}.meta; end

merged.meta = hcat_structs(metas{:});











function merged = hcat_structs(varargin)

% Recursively figure out all the fields present in each struct.
ds_fields = repmat({{}}, 1, length(varargin));
for d = 1:length(varargin)
	ds = varargin{d};
	fstack = {''};
	while ~isempty(fstack)
		fields = eval(sprintf('fieldnames(ds%s)', fstack{1}));
		for f = 1:length(fields)
			nf = [fstack{1} '.' fields{f}];
			eval(sprintf('is_struct = isstruct(ds%s);', nf));
			if is_struct
				fstack{end+1} = nf;
			else
				ds_fields{d}{end+1} = nf;
			end
		end
		fstack = fstack(2:end);
	end
end

fields = unique(cat(2, ds_fields{:}));
fields = fields(~rx(fields, '^\.type$'));




% Calculate the number of samples in each dataset.
S = nan(1, length(varargin));
for d = 1:length(varargin)
	S(d) = length(varargin{d}.sample_id);
end
S_ranges = [[1; cumsum(S(1:end-1))'+1], [cumsum(S)']];
S_total = sum(S);



% Construct a FxD matrix where each row indicates which datasets have the field
% and what the field's type is (0 = missing, 1 = numeric, 2 = cell).
field_type = zeros(length(fields), length(varargin));
for d = 1:length(varargin)
	field_type(:, d) = ismember(fields, ds_fields{d});
	for k = find(field_type(:, d))'
		field_type(k, d) = 1+eval(sprintf('iscell(varargin{d}%s)', fields{k}));
	end
end





% Merge the datasets field by field.
merged = struct;
for f = 1:length(fields)
	if any(field_type(f, :) == 2)
		eval(sprintf('merged%s = repmat({''''}, 1, S_total);', fields{f}));
		for d = find(field_type(f, :) > 0)
			eval(sprintf( ...
				'merged%s(S_ranges(d,1):S_ranges(d,2)) = varargin{d}%s;', ...
				fields{f}, fields{f}));
		end
	else
		eval(sprintf('merged%s = nan(1, S_total);', fields{f}));
		for d = find(field_type(f, :) > 0)
			eval(sprintf( ...
				'merged%s(S_ranges(d,1):S_ranges(d,2)) = varargin{d}%s;', ...
				fields{f}, fields{f}));
		end
	end
end







function [jlabel, jmatrix] = join_matrices(labels, matrices)

% Double check that the labels are unique in each matrix.
for d = 1:length(matrices)
	if length(unique(labels{d})) ~= length(labels{d})
		error('Labels not unique in dataset #%d.', d);
	end
end

all_labels = unique(cat(1, labels{:}));
for d = 1:length(matrices)
	[~, label_pos{d}] = ismember(all_labels, labels{d});
end

% Count the total number of columns.
S = 0;
for d = 1:length(matrices), S = S + size(matrices{d}, 2); end

jlabel = all_labels;
jmatrix = nan(length(all_labels), S);

s = 0;
for d = 1:length(matrices)
	jmatrix(label_pos{d} > 0, s+1:s+size(matrices{d}, 2)) = ...
		matrices{d}(label_pos{d}(label_pos{d} > 0), :);
	s = s + size(matrices{d}, 2);
end

