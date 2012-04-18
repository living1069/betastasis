function [] = gset_diff_expressed(A, B, gset, varargin)

global organism;

if ~isempty(B) && size(A, 1) ~= size(B, 1)
	error 'Matrices for classes A and B must have an equal number of rows.';
end

valid = find(~any(isnan(A), 2) & ~any(isnan(B), 2));
gset = intersect(gset, valid);

diff = log2(median(A(gset, :), 2) ./ median(B(gset, :), 2));

[~, order] = sort(diff, 'descend');

if length(order) > 50
	order = order([1:25, length(order)-24:length(order)]);
end

diff = diff(order);
gset = gset(order);

diff = diff(:);
gene_labels = organism.Genes.Name(gset);

diff = [diff(1:25); 0; diff(26:50)];
gene_labels = cat(1, gene_labels(1:25), {'...'}, gene_labels(26:50));

%w = 50;

figure('PaperUnits', 'normalized', 'PaperPosition', [-0.05 0.2 1.1 0.7]);
shading flat; bar(diff, 'black'); ylabel('Logratio'); xlim([0 52]);
set(gca, 'xtick', 1:length(gene_labels), 'xticklabel', gene_labels);
set(gca, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
rotateticklabel(gca, 90);
saveas(gcf, '~/gset_diff.pdf');




