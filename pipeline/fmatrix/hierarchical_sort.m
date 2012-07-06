function order = hierarchical_sort(varargin)

sort_score = zeros(1, length(varargin{1}));

for k = 1:length(varargin)
	tmp = varargin{k} - min(varargin{k}) + 1;
	tmp(isnan(tmp)) = 0;
	sort_score = sort_score * max(tmp);
	sort_score = sort_score + tmp;
end

[~, order] = sort(sort_score, 'descend');

