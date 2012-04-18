
% Author: Matti Annala <matti.annala@tut.fi>

function [] = find_fusions_exon_array(exon_expr, varargin)
	
global organism;
genes = organism.Genes;
exons = organism.Exons;

sample_ids = exon_expr.Meta.Sample.ID;

S = size(exon_expr.Mean, 2);

t_stats = nan(length(genes.Name), S);

for g = 1:length(genes.Name)
	gene_exons = find(exons.Gene == g);

	[~, order] = sort_nat(exons.ID(gene_exons));
	gene_exons = gene_exons(order);
	if isempty(gene_exons), continue, end

	E = length(gene_exons);

	gex_expr = exon_expr.Mean(gene_exons, :) + 1;
	gene_expr = nansum(gex_expr, 1) + 2;
	
	if isnan(gene_expr), continue, end

	test_ratio = gex_expr ./ repmat(gene_expr, E, 1);
	test_ratio = test_ratio ./ (1 - test_ratio);
	test_ratio = log2(test_ratio);
	
	test_ratio = test_ratio - repmat(nanmedian(test_ratio, 2), 1, S);
	
	% Here we assume that the odds ratios are roughly normally distributed.
	stdev = nanstd(test_ratio(:));
	test_ratio = test_ratio / stdev;
	
	if g == gene_idx('TACC3')
		figure; ksdensity(test_ratio(:));
		saveas(gcf, 'TACC3_odds_ratios.pdf');
	else
		continue;
	end
	
	% Calculate the T test statistic for each sample specifically.
	S_r = zeros(1, S);
	S_n = nansum(test_ratio, 1);
	for r = 1:E-1
		S_r = S_r + test_ratio(r, :);
		t_stats(g, :) = max(t_stats(g, :), ...
			(S_r / r - (S_n - S_r) ./ (E - r)) ./ sqrt(1/r + 1/(E-r)));
	end
end

t_stats_ns = t_stats(:);
valid = find(~isnan(t_stats_ns));

[~, order] = sort(t_stats_ns(valid), 'descend');
[g, s] = ind2sub(size(t_stats), valid(order));

whos order g s

for k = 1:100
	fprintf(1, 'Gene %s in sample %s.\n', genes.Name{g(k)}, sample_ids{s(k)});
end







