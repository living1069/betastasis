function [row_perm, col_perm] = better_cluster(data, varargin)

dist_metric = 'euclidean';
linkage_method = 'average';
clim = [];

if ~isempty(varargin)
	for k = 1:2:length(varargin)
		if strcmpi(varargin{k}, 'CLim')
			clim = varargin{k+1};
			continue;
		end
		
		if strcmpi(varargin{k}, 'DistanceMetric')
			dist_metric = varargin{k+1};
			continue;
		end
		
		if strcmpi(varargin{k}, 'Linkage')
			linkage_method = varargin{k+1};
			continue;
		end

		error('Unrecognized option "%s".', varargin{k});
	end
end

% First we cluster by columns.
figure;
clust = linkage(pdist(data', dist_metric), linkage_method);
[~, ~, col_perm] = dendrogram(clust, 0, 'ColorThreshold', 'default');
data = data(:, col_perm);
saveas(gcf, '~/cluster_column_dendro.pdf');

% Then we cluster by rows.
figure;
clust = linkage(pdist(data, dist_metric), linkage_method);
[~, ~, row_perm] = dendrogram(clust, 0, 'ColorThreshold', 'default');
data = data(row_perm, :);
saveas(gcf, '~/cluster_row_dendro.pdf');







% Render the clustered heatmap.

%[~, order] = sort(data(:));
%rank = 1:length(order);
%rank(order) = rank;
%data(:) = rank;

%range = max(-min(min(data)), max(max(data)))
%whos data
%data = round((data + range) / (2*range) * 255);
%expand = 2;

%min(min(data))
%max(max(data))

%imwrite(data / numel(data) * 255, jet, '~/cluster.png', 'png');

figure; colormap(redgreencmap(256));
if isempty(clim)
	imagesc(data);
else
	imagesc(data, clim);
end
set(gca, 'Visible', 'off');
saveas(gcf, '~/cluster_heatmap.png');









function [] = binarize()

for g = 1:size(b, 1)
	binary = kmeans(b(g, :)', 2, 'Start', [min(b(g, :)); max(b(g, :))]);
	b(g, :) = binary' - 1;
end

