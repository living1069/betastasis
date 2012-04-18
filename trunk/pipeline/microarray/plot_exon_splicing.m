
% Author: Matti Annala <matti.annala@tut.fi>

function [] = plot_exon_splicing(gene, test_exon, ref_exon)

global organism;
genes = organism.Genes;
exons = organism.Exons;

gene = gene_idx(gene);
gene_exons = find(exons.Gene == gene);

[~, order] = sort_nat(exons.ID(gene_exons));
gene_exons = gene_exons(order);

St = size(test_exon.Mean, 2);
Sr = size(ref_exon.Mean, 2);
E = length(gene_exons);

test_exon_expr = test_exon.Mean(gene_exons, :) + 1;
ref_exon_expr = ref_exon.Mean(gene_exons, :) + 1;
test_gene_expr = nansum(test_exon_expr, 1) + 2;
ref_gene_expr = nansum(ref_exon_expr, 1) + 2;

if 1
	% Normalize exon expressions based on the expression of other exons of the
	% same gene. The equations become analogous to odds ratio calculations
	% in statistics.
	test_ratio = test_exon_expr ./ repmat(test_gene_expr, E, 1);
	test_ratio = test_ratio ./ (1 - test_ratio);
	test_ratio = log2(test_ratio);

	ref_ratio = ref_exon_expr ./ repmat(ref_gene_expr, E, 1);
	ref_ratio = ref_ratio ./ (1 - ref_ratio);
	ref_ratio = log2(ref_ratio);
	
	test_ratio = test_ratio - repmat(nanmedian(ref_ratio, 2), 1, St);
	ref_ratio = ref_ratio - repmat(nanmedian(ref_ratio, 2), 1, Sr);
else
	% Just normalize based on median exon expression levels.
	test_ratio = test_exon_expr ./ repmat(nanmedian(ref_exon_expr, 2), 1, St);
	test_ratio = log2(test_ratio);
end

valid = find(any(isnan(test_ratio), 2));

figure; hold all;

for s = 1:Sr
	plot(valid, ref_ratio(valid, s), 'Color', [.9 .9 .9]);
end

colors = [1 0 0; 0 1 0; 0 0 1];
for s = 1:St
	plot(valid, test_ratio(valid, s), 'Color', colors(s, :));
end

set(gca, 'XTick', valid, 'XTickLabel', exons.ID(gene_exons(valid)));
xlabel('Exon ID'); ylabel('Log odds ratio');
saveas(gcf, 'exon_splice.pdf');


