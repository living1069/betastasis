function [order, sort_score] = hierarchical_sort(varargin)

direction = 'ascend';

descend = strcmpi('descend', varargin);
if any(descend)
	direction = 'descend';
	varargin = varargin(~descend);
end

sort_score = zeros(1, length(varargin{1}));

for k = 1:length(varargin)
	tmp = varargin{k} - min(varargin{k}) + 2;
	tmp = tmp(:)';
	tmp(isnan(tmp)) = 1;
	sort_score = sort_score * max(tmp);
	sort_score = sort_score + tmp;
end

[~, order] = sort(sort_score, direction);

