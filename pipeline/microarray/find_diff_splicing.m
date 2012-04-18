
% FIND_DIFF_SPLICING    Find differential splicing events from expression data

% Author: Matti Annala <matti.annala@tut.fi>

function [] = find_diff_splicing(test_exon, test_gene, ref_exon, ref_gene, varargin)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;

test_exon_log_expr = log2(test_exon.Mean);
test_gene_log_expr = log2(test_gene.Mean);
ref_exon_log_expr = log2(ref_exon.Mean);
ref_gene_log_expr = log2(ref_gene.Mean);

if size(test_exon_log_expr, 2) ~= size(test_gene_log_expr, 2) || ...
	size(ref_exon_log_expr, 2) ~= size(ref_gene_log_expr, 2)
	error 'Exon and gene expression matrices must have equal dimensions.';
end

if ~all(strcmpi(test_exon.Meta.Sample.ID, test_gene.Meta.Sample.ID))
	fprintf(1, 'WARNING: Sample IDs do not match for the test samples.\n');
end
if ~all(strcmpi(ref_exon.Meta.Sample.ID, ref_gene.Meta.Sample.ID))
	fprintf(1, 'WARNING: Sample IDs do not match for the reference samples.\n');
end

test_splice_index = nan(size(test_exon_log_expr));
ref_splice_index = nan(size(ref_exon_log_expr));

E = size(test_exon_log_expr, 1);

for k = 1:E
	g = exons.Gene(k);
	test_splice_index(k, :) = ...
		test_exon_log_expr(k, :) - test_gene_log_expr(g, :);
	ref_splice_index(k, :) = ...
		ref_exon_log_expr(k, :) - ref_gene_log_expr(g, :);
end

test_sindex_mean = nanmedian(test_splice_index, 2);
ref_sindex_mean = nanmedian(ref_splice_index, 2);

[~, p] = ttest2(test_splice_index', ref_splice_index');
aberrated = nansum(abs(test_splice_index - ...
	repmat(ref_sindex_mean, 1, size(test_splice_index, 2))) > 3.0, 2);
valid = find(~isnan(aberrated));
[~, order] = sort(aberrated(valid), 'descend');
order = valid(order);

fprintf(1, 'Gene\tExon\tAberrated (%%)\tP-value\tDiff. splice\tGene logratio\n');

order = order(1:100);


diff_genes = false(length(genes.Name), 1);
for k = 1:length(order)
	ex = order(k);
	g = exons.Gene(ex);
	
	if diff_genes(g), continue, end
	diff_genes(g) = true;
	
	gene_exons = find(exons.Gene == g);
	diff_gene_exons = intersect(gene_exons, order);
	
	[~, nat_order] = sort_nat(exons.ID(diff_gene_exons));
	diff_gene_exons = diff_gene_exons(nat_order);
	
	for x = 1:length(diff_gene_exons)
		ex = diff_gene_exons(x);
		fprintf(1, '%s\t%s\t%.1f%%\t%g\t%.2f\t%.2f\n', ...
			genes.Name{g}, exons.ID{ex}, ...
			aberrated(ex) / size(test_splice_index, 2) * 100, p(ex), ...
			test_sindex_mean(ex) - ref_sindex_mean(ex), ...
			nanmedian(test_gene_log_expr(g, :), 2) - ...
			nanmedian(ref_gene_log_expr(g, :), 2));
	end
end


