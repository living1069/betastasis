function [row_perm, col_perm] = better_cluster(data, varargin)

cluster_rows = true;
cluster_columns = true;
dist_metric = 'euclidean';
linkage_method = 'average';

keep = true(1, length(varargin));
for k = 1:2:length(varargin)
	if rx(varargin{k}, 'metric')
		dist_metric = varargin{k+1}; keep(k:k+1) = false;
	end
	if rx(varargin{k}, 'link')
		linkage_method = varargin{k+1}; keep(k:k+1) = false;
	end
	if rx(varargin{k}, 'clust.*rows')
		cluster_rows = varargin{k+1}; keep(k:k+1) = false;
	end
	if rx(varargin{k}, 'clust.*col')
		cluster_columns = varargin{k+1}; keep(k:k+1) = false;
	end
end
varargin = varargin(keep);

% First we cluster by columns.
if cluster_columns
	figure;
	clust = linkage(pdist(data', dist_metric), linkage_method);
	col_perm = dendroperm(clust, 0);
	data = data(:, col_perm);
	saveas(gcf, '~/cluster_column_dendro.pdf');
else
	col_perm = 1:size(data, 2);
end

% Then we cluster by rows.
if cluster_rows
	figure;
	clust = linkage(pdist(data, dist_metric), linkage_method);
	row_perm = dendroperm(clust, 0);
	data = data(row_perm, :);
	saveas(gcf, '~/cluster_row_dendro.pdf');
else
	row_perm = 1:size(data, 1);
end

heatmap(data, '~/cluster_heatmap.png', varargin{:});

