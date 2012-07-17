
% PAIR_ROWS    Associate the rows of two or more datasets
%
%    This function works the same way as the function PAIR, with the exception
%    that this functions acts on dataset rows. See HELP PAIR for documentation.

% Author: Matti Annala <matti.annala@tut.fi>

function varargout = pair_rows(varargin)

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
	keys{k} = eval(['varargin{k}.rows.' key ';']);
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
		varargout{k} = filter_rows(varargin{k}, pos);
	end
	
elseif rx(method, 'union')
	uniq = unique(vertcat(keys{:}));
	
	for k = 1:length(varargin)-1
		[~, pos] = ismember(uniq, keys{k});
		pos(pos == 0) = NaN;
		varargout{k} = nanfilter_rows(varargin{k}, pos);
	end
end


