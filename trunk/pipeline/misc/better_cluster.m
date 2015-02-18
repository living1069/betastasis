function [row_perm, col_perm] = better_cluster(data, varargin)

cluster_rows = true;
cluster_columns = true;
dendro_rows = false;
dendro_columns = false;
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
	if rx(varargin{k}, 'dendro.*row')
		dendro_rows = varargin{k+1}; keep(k:k+1) = false;
	end
	if rx(varargin{k}, 'dendro.*col')
		dendro_columns = varargin{k+1}; keep(k:k+1) = false;
	end
end
varargin = varargin(keep);

if cluster_columns
	fprintf('Clustering columns...\n');
	clust = linkage(pdist(data', dist_metric), linkage_method);
	if dendro_columns
		figure; [~, ~, col_perm] = dendrogram(clust, 0);
		saveas(gcf, '~/dendro_column.pdf');
	else
		col_perm = dendroperm(clust, 0);
	end
	data = data(:, col_perm);
else
	col_perm = 1:size(data, 2);
end

if cluster_rows
	fprintf('Clustering rows...\n');
	clust = linkage(pdist(data, dist_metric), linkage_method);
	if dendro_rows
		figure; [~, ~, row_perm] = dendrogram(clust, 0);
		saveas(gcf, '~/dendro_rows.pdf');
	else
		row_perm = dendroperm(clust, 0);
	end
	data = data(row_perm, :);
else
	row_perm = 1:size(data, 1);
end

heatmap(data, '~/cluster_heatmap.png', varargin{:});

