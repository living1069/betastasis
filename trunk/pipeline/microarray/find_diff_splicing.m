
% FIND_DIFF_SPLICING    Find differential splicing events from expression data

% Author: Matti Annala <matti.annala@tut.fi>

function [] = find_diff_splicing(test_exon, test_gene, ref_exon, ref_gene, ...
	varargin)

global organism;
genes = organism.Genes;
exons = organism.Exons;


test_exon_log_expr = log2(test_exon.Mean);
test_gene_log_expr = log2(test_gene.Mean);
ref_exon_log_expr = log2(ref_exon.Mean);
ref_gene_log_expr = log2(ref_gene.Mean);

if size(test_exon_log_expr, 2) ~= size(test_gene_log_expr, 2) || ...
	size(ref_exon_log_expr, 2) ~= size(ref_gene_log_expr, 2)
	error 'Exon and gene expression matrices must have equal dimensions.';
end

test_splice_index = nan(size(test_exon_log_expr));
ref_splice_index = nan(size(ref_exon_log_expr));

E = size(test_exon_log_expr, 1);
significance = 0.05 / E

for k = 1:E
	g = exons.Gene(k);
	test_splice_index(k, :) = ...
		test_exon_log_expr(k, :) - test_gene_log_expr(g, :);
	ref_splice_index(k, :) = ...
		ref_exon_log_expr(k, :) - ref_gene_log_expr(g, :);
end

[~, p] = ttest2(test_splice_index', ref_splice_index');

test_sindex_mean = nanmean(test_splice_index, 2);
ref_sindex_mean = nanmean(ref_splice_index, 2);

sig = find(p <= significance);

[~, order] = sort(p(sig), 'ascend');
sig = sig(order);

fprintf(1, 'Differentially spliced exons:\n');
for k = 1:length(sig)
	ex = sig(k);
	fprintf(1, '- %s[%s], p-value %g (%.2f vs %.2f)\n', ...
		genes.Name{exons.Gene(ex)}, exons.ID{ex}, p(ex), ...
		test_sindex_mean(ex), ref_sindex_mean(ex));
end


