
% Author: Matti Annala <matti.annala@tut.fi>

function [] = find_diff_splicing_rnaseq(test_exon, ref_exon, varargin)

global organism;
genes = organism.Genes;
exons = organism.Exons;

S = size(test_exon.Mean, 2);

test_gene_expr = nan(length(genes.Name), S);
ref_gene_expr = nan(length(genes.Name), size(ref_exon.Mean, 2));

for g = 1:length(genes.Name)
	gex = find(exons.Gene == g);
	test_gene_expr(g, :) = sum(test_exon.Mean(gex, :), 1);
	ref_gene_expr(g, :) = sum(ref_exon.Mean(gex, :), 1);
end

test_exon_expr = test_exon.Mean + 1;
ref_exon_expr = ref_exon.Mean + 1;
test_gene_expr = test_gene_expr + 2;
ref_gene_expr = ref_gene_expr + 2;

E = size(test_exon_expr, 1);

test_odds_ratio = test_exon_expr ./ test_gene_expr(exons.Gene, :);
test_odds_ratio = test_odds_ratio ./ (1 - test_odds_ratio);
test_odds_ratio = log2(test_odds_ratio);

ref_odds_ratio = ref_exon_expr ./ ref_gene_expr(exons.Gene, :);
ref_odds_ratio = ref_odds_ratio ./ (1 - ref_odds_ratio);
ref_odds_ratio = log2(ref_odds_ratio);

test_sindex_mean = median(test_odds_ratio, 2);
ref_sindex_mean = median(ref_odds_ratio, 2);

[~, p] = ttest2(test_odds_ratio', ref_odds_ratio');
aberrated = nansum(abs(test_odds_ratio - ...
	repmat(ref_sindex_mean, 1, size(test_odds_ratio, 2))) > 4.0, 2);

ranking = 'copa';
if regexpi(ranking, 'aberrated')
	[~, order] = sort(aberrated, 'descend');
elseif regexpi(ranking, 'ttest')
	[~, order] = sort(p, 'ascend');
elseif regexpi(ranking, 'copa')
	mad = median(abs(ref_odds_ratio - ...
		repmat(ref_sindex_mean, 1, size(ref_odds_ratio, 2))), 2);
	valid = find(mad > 0);
	
	score = quantile(abs(test_odds_ratio - ...
		repmat(ref_sindex_mean, 1, size(test_odds_ratio, 2)) - ...
		repmat(mad, 1, size(test_odds_ratio, 2))), 0.75, 2);
	
	[~, order] = sort(score(valid), 'descend');
	order = valid(order);
	
	score(order(1:10))
end

fid = fopen('top_diff_splicing_rnaseq.txt', 'W');
fprintf(fid, 'Gene\tExon\tAberrated (%%)\tP-value\tLog-OR\tGene logratio\n');

order = order(1:1000);

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
		fprintf(fid, '%s\t%s\t%.1f%%\t%g\t%.2f\t%.2f\n', ...
			genes.Name{g}, exons.ID{ex}, ...
			aberrated(ex) / size(test_odds_ratio, 2) * 100, p(ex), ...
			test_sindex_mean(ex) - ref_sindex_mean(ex), ...
			log2(nanmedian(test_gene_expr(g, :), 2) / ...
			nanmedian(ref_gene_expr(g, :), 2)));
	end
end

fclose(fid);

