function [order, sort_score] = hierarchical_sort(varargin)

sort_score = zeros(1, length(varargin{1}));

for k = 1:length(varargin)
	tmp = varargin{k} - min(varargin{k}) + 2;
	tmp = tmp(:)';
	tmp(isnan(tmp)) = 1;
	sort_score = sort_score * max(tmp);
	sort_score = sort_score + tmp;
end

[~, order] = sort(sort_score);

