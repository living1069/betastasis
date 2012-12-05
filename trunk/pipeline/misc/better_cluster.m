function [row_perm, col_perm] = better_cluster(data, varargin)

cluster_rows = true;
dist_metric = 'euclidean';
linkage_method = 'average';
color_map = redgreencmap(256);
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
		
		if rx(varargin{k}, 'color')
			color_map = varargin{k+1};
			continue;
		end
		
		if rx(varargin{k}, 'cluster.*rows')
			cluster_rows = varargin{k+1};
			continue;
		end

		error('Unrecognized option "%s".', varargin{k});
	end
end

% First we cluster by columns.
figure;
d = pdist(data', dist_metric);
clust = linkage(pdist(data', dist_metric), linkage_method);
[~, ~, col_perm] = dendrogram(clust, 0, 'ColorThreshold', 'default');
data = data(:, col_perm);
saveas(gcf, '~/cluster_column_dendro.pdf');

% Then we cluster by rows.
if cluster_rows
	figure;
	clust = linkage(pdist(data, dist_metric), linkage_method);
	[~, ~, row_perm] = dendrogram(clust, 0, 'ColorThreshold', 'default');
	data = data(row_perm, :);
	saveas(gcf, '~/cluster_row_dendro.pdf');
else
	row_perm = 1:size(data, 1);
end




im = data;
if isempty(clim), clim = [nanmin(data(:)), nanmax(data(:))]; end
if length(clim) == 1, clim = [-clim 0 clim];
elseif length(clim) == 2, clim = [clim(1) (clim(1)+clim(2))/2 clim(2)];
end


if size(color_map, 1) == 256
	fprintf('Colormap contains 256 colors. Interpolating data values...');
	im(:) = interp1(clim, [1 128 256], im(:), 'linear', 'extrap');
	im(im < 1) = 1;
	im(im > 256) = 256;
elseif all(all(round(im) == im | isnan(im))) && nanmin(im(:)) > 0
	fprintf('Data values are indices. Indexing into specified color map...');
	color_map(end+1, :) = [1 1 1]; nan_idx = length(color_map);
	im(isnan(im)) = nan_idx;
else
	error('ERROR: Not sure how to render colors.');
end

imwrite(im, color_map, '~/cluster_heatmap.png');


