function labels = cluster_expr(expr, varargin)

dist_metric = 'euclidean';
linkage_method = 'average';

if ~isempty(varargin)
	for k = 1:2:length(varargin)
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

log_expr = log2(expr.Mean);

% Prefilter out genes whose expression levels are unknown or that have a very
% low standard deviation.
not_available = (sum(isnan(log_expr), 2) > 0);
low_sdev = (std(log_expr, 0, 2) < 1);

included_genes = find(~(not_available | low_sdev));
length(included_genes)

c = log_expr(included_genes, :);

b = c;
for g = 1:size(b, 1)
	binary = kmeans(b(g, :)', 2, 'Start', [min(b(g, :)); max(b(g, :))]);
	b(g, :) = binary' - 1;
end

%dist = squareform(pdist(b', 'cityblock'));
dist = squareform(pdist(c', 'euclidean'));
sep_score = zeros(size(b, 1), 1);

for g = 1:size(b, 1)
	binary = b(g, :);
	
	low = (binary == 0);
	high = (binary == 1);
	
	if min(sum(low), sum(high)) < 5
		sep_score(g) = -Inf;
	else
		sep_score(g) = mean(mean(dist(low, high)));
	end
end

mu = mean(sep_score(sep_score ~= -Inf));
sigma = std(sep_score(sep_score ~= -Inf));
sep_p_score = normcdf(sep_score, mu, sigma);

best = (sep_p_score > 0.99);
%sum(best)

worst = (sep_p_score < 0.50);
sum(~worst)

c = c(~worst, :);



dist = pdist(c', 'euclidean');
clust = linkage(dist, 'average');

figure;
[~, ~, perm] = dendrogram(clust, 0, 'ColorThreshold', 'default');
saveas(gcf, '~/dendro.pdf');

labels = cluster(clust, 'Cutoff', 0.7 * max(clust(:,3)), ...
	'Criterion', 'distance');
labels = labels(perm);

% Permute the labels so that they agree with the dendrogram leaf order.
group_count = 0;
new_labels = zeros(size(labels));
for k = 1:length(labels)
	if new_labels(k) == 0
		group_count = group_count + 1;
		new_labels(labels == labels(k)) = group_count;
	end
end

labels(perm) = new_labels;

for k = 1:max(labels)
	fprintf(1, 'Cluster %d: %d items\n', k, sum(labels == k));
end



p = c(:, perm);

dist = pdist(p, 'euclidean');
clust = linkage(dist, 'average');
[~, ~, perm] = dendrogram(clust, 0);

p = p(perm, :);

%bandh = round(size(p, 1) / 30);
%band_bg_color = 0.5*min(min(p)) + 0.5*max(max(p));
%band_fg_color = max(max(p));
%p = [band_bg_color * ones(bandh, size(p, 2)); p];

figure; colormap(redgreencmap(256));
imagesc(p); saveas(gcf, '~/cluster.png');


