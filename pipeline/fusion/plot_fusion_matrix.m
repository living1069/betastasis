function [] = plot_fusion_matrix(fusion_genes, fusions)

global organism;
genes = organism.Genes;

dist_metric = 'cityblock';
linkage_method = 'average';

S = length(fusions.Fusions);

if ~iscellstr(fusion_genes)
	error 'Fusion genes must be specified as a cell array of strings.';
end

fusion_matrix = zeros(length(fusion_genes), S);
for k = 1:length(fusion_genes)
	fusion_matrix(k, :) = reads_for_fusion(fusion_genes{k}, fusions) > 0;
end



if 0 
figure;
clust = linkage(pdist(fusion_matrix', dist_metric), linkage_method);
[~, ~, col_perm] = dendrogram(clust, 0, 'ColorThreshold', 'default');
data = fusion_matrix(:, col_perm);
saveas(gcf, '~/cluster_column_dendro.pdf');

figure; colormap(redgreencmap(256));
imagesc(data); set(gca, 'Visible', 'off');
saveas(gcf, '~/cluster_heatmap.png');

else

fixed = 0;
for f = 1:length(fusion_genes)
	pivot = fusion_matrix(:, fixed+1:end);
	new_fixed = sum(pivot(f, :));
	fusion_matrix(:, fixed+1:end) = ...
		pivot(:, [find(pivot(f, :)) find(~pivot(f, :))]);
	fixed = fixed + new_fixed;
end

figure; colormap(redgreencmap(256));
imagesc(fusion_matrix); set(gca, 'Visible', 'off');
text(repmat(S-1, 1, length(fusion_genes)), 1:length(fusion_genes), ...
	fusion_genes, 'HorizontalAlignment', 'right');
saveas(gcf, '~/cluster_heatmap.png');

end

