function [labels, clustered_order] = cluster_cgh(samples, refs, ...
	cgh_probesets, varargin)

save_dendrogram = '';
region_chr = 0;
region_pos = [];
distance_metric = 'euclidean';

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Dendrogram')
		save_dendrogram = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
	
	if strcmpi(varargin{k}, 'Region')
		tokens = regexpi(varargin{k+1}, '(chr.+?):\s*(\d+)\s*-\s*(\d+)', ...
			'tokens');
		if length(tokens) ~= 1, error 'Invalid region specification.'; end
		tokens = tokens{1};
		region_chr = chromosome_sym2num(tokens{1});
		region_pos = [str2double(tokens{2}) str2double(tokens{3})];
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

cna = cna_from_cgh(samples, refs, cgh_probesets, varargin{:});

if region_chr
	region_ps = find(cgh_probesets.Chromosome == region_chr & ...
		cgh_probesets.Offset >= region_pos(1) & ...
		cgh_probesets.Offset <= region_pos(2));
	if any(region_ps ~= (region_ps(1):region_ps(end))')
		error 'Probesets are not ordered by chromosome.';
	end
	
	cna = cna(region_ps, :);
end

bad_loci = any(isnan(cna), 2);
cna = cna(~bad_loci, :);

dist = pdist(cna', distance_metric);
clust = linkage(dist, 'average');

figure; [~, ~, perm] = dendrogram(clust, 0, 'ColorThreshold', 'default');
if ~isempty(save_dendrogram)
	saveas(gcf, save_dendrogram);
end

clustered_order = perm;

labels = cluster(clust, 'Cutoff', 0.7 * max(clust(:, 3)), ...
	'Criterion', 'distance');

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






