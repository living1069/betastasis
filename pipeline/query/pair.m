
% PAIR      Associate the columns of two or more datasets
%
%    [PA, PB] = PAIR(A, B, KEY) prepares two datasets A and B for integrative
%    analysis, returning two new datasets PA and PB with an equal number of
%    columns. The columns in PA and PB are ordered according to the meta
%    field KEY, and any columns not shared by both datasets A and B are
%    discarded.
%
%    This function is useful in experiments where multiple test samples need to
%    be compared with their respective reference samples in a pairwise fashion.
%
%    Example:
%    We have a dataset A containing three tumor samples from patients
%    A23, A54 and A12. We also have a dataset B containing two adjacent normal
%    samples from patients A12 and A54. We run
%    
%        [PA, PB] = PAIR(A, B, 'patient_id')
%    
%    to build two new datasets PA and PB where the samples are paired
%    according to patient ID. The columns in datasets PA and PB now correspond
%    to patients A12 and A54 (note that the original order was lost).
%    The tumor sample for patient A23 was discarded, since no
%    matching reference sample was found in dataset B.

% Author: Matti Annala <matti.annala@tut.fi>

function varargout = pair(varargin)

method = 'intersect';
if ischar(varargin{end}) && rx(varargin{end}, 'union')
	method = 'union';
	varargin = varargin(1:end-1);
end

if ~ischar(varargin{end})
	error 'The last argument must be a string that specifies the join key.';
end

key = varargin{end};
for k = 1:length(varargin)-1
	keys{k} = eval(['varargin{k}.meta.' key ';']);
	if length(unique(keys{k})) ~= length(keys{k})
		error('Key in dataset #%d is not one-to-one.', k);
	end
end

if rx(method, 'intersect')
	uniq = intersect(keys{1}, keys{2});
	for k = 3:length(keys)
		uniq = intersect(uniq, keys{k});
	end

	for k = 1:length(varargin)-1
		[~, pos] = ismember(uniq, keys{k});
		varargout{k} = filter(varargin{k}, pos);
	end
	
elseif rx(method, 'union')
	uniq = unique(horzcat(keys{:}));
	
	for k = 1:length(varargin)-1
		[~, pos] = ismember(uniq, keys{k});
		pos(pos == 0) = NaN;
		varargout{k} = nanfilter(varargin{k}, pos);
	end
end

