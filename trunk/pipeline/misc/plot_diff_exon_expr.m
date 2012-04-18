function [] = plot_diff_exon_expr(test, ref, gene)

global organism;
genes = organism.Genes;
exons = organism.Exons;

if ischar(gene), gene = gene_idx(gene); end

gene_exons = find(exons.Gene == gene);

test_gene_expr = nansum(test.Mean(gene_exons, :), 1);
ref_gene_expr = nansum(ref.Mean(gene_exons, :), 1);

test_exon_expr = log2(test.Mean(gene_exons, :)) - ...
	repmat(log2(test_gene_expr), length(gene_exons), 1);
ref_exon_expr = log2(ref.Mean(gene_exons, :)) - ...
	repmat(log2(ref_gene_expr), length(gene_exons), 1);

ref_exon_expr = median(ref_exon_expr, 2);

expr_logratio = test_exon_expr - ...
	repmat(ref_exon_expr, 1, size(test_exon_expr, 2));
	
for s = 1:size(expr_logratio, 2)
	figure; plot(expr_logratio(:, s));
	ylim([-5 5]);
	saveas(gcf, sprintf('%s.pdf', test.Meta.Sample.ID{s}));
end


