
% META_TABULAR    Augment dataset with metadata from a TLM or XLS file.
%    
%    META = META_TABULAR(TLM, META) reads patient and sample specific metadata
%    from the tab delimited file TLM, and associates the metadata with
%    corresponding samples in query set META.

% Author: Matti Annala <matti.annala@tut.fi>

function meta = meta_tabular(filepath, meta, varargin)

join_field = 'sample_id';

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'join')
		join_field = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end


if rx(filepath, '\.(txt|csv|tsv)$')
	[data, headers] = readtable(filepath);
	
	% Remove columns without a proper header.
	bad_cols = find(cellfun(@isempty, headers));
	if ~isempty(bad_cols), bad_cols, end
	headers(bad_cols) = [];
	raw(:, bad_cols) = [];
	
	% If a column contains nothing but numbers and empty cells, convert it
	% into a numeric vector.
	for k = 1:length(data)
		empty = strcmp(data{k}, '') | strcmpi(data{k}, 'NaN');
		num = str2double(data{k}(~empty));
		if all(~isnan(num))
			data{k} = str2double(data{k});
		end
	end

else
	error 'Unsupported metadata format.';
end

% First we build a metadata structure with all of the data.
new_meta = struct;
new_meta_fields = {};
for k = 1:length(headers)
	fname = lower(headers{k});
	fname = regexprep(fname, '[ -]', '_');
	fname = regexprep(fname, '[():]', '');
	fname = regexprep(fname, '[_.]+$', '');
	fname = strrep(fname, '%', 'percent');
	new_meta_fields{k} = fname;
end

sample_id_fields = new_meta_fields(rx(new_meta_fields, '^sample.*(id|name)'));

if length(sample_id_fields) < 1
	error 'New metadata does not contain a sample ID column.';
elseif length(sample_id_fields) > 1
	error 'Too many candidate sample ID columns found.';
end

for k = 1:length(new_meta_fields)
	if strcmp(new_meta_fields{k}, sample_id_fields{1})
		new_meta_fields{k} = 'sample_id';
		new_meta.sample_id = data{k}';
	else
		eval(sprintf('new_meta.%s = data{k}'';', new_meta_fields{k}));
	end
end

if length(unique(new_meta.sample_id)) ~= length(new_meta.sample_id)
	new_meta.sample_id'
	error 'New metadata must not contain duplicate sample IDs.';
end

% And that's it...
if nargin == 1
	meta = new_meta;
	return;
end

% ... unless the user wants to associate the new metadata with an existing
% dataset. In that case, construct a mapping between the old and new metadata,
% based on the sample ID.
if ~isfield(meta, 'sample_id')
	error 'Old metadata does not contain a sample ID.';
end

S = length(meta.sample_id);
%eval(['meta.' join_field])'
%eval(['new_meta.' join_field])'
[~, pos] = ismember(eval(['meta.' join_field]), eval(['new_meta.' join_field]));

for k = 1:length(new_meta_fields)
	if strcmp(new_meta_fields{k}, join_field), continue, end
		
	nf = eval(sprintf('new_meta.%s', new_meta_fields{k}));
	
	if iscellstr(nf)
		eval(sprintf('meta.%s = repmat({''''}, 1, S);', new_meta_fields{k}));
		eval(sprintf('meta.%s(pos ~= 0) = nf(pos(pos ~= 0));', ...
			new_meta_fields{k}));
	elseif isnumeric(nf)
		eval(sprintf('meta.%s = nan(1, S);', new_meta_fields{k}));
		eval(sprintf('meta.%s(pos ~= 0) = nf(pos(pos ~= 0));', ...
			new_meta_fields{k}));
	else
		error('Field %s is of unsupported type.', new_meta_fields{k});
	end
end

