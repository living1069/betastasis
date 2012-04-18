
% Author: Matti Annala <matti.annala@tut.fi>

function enriched = hypergeotest(white, gsets)

global organism;

if ~islogical(white), error 'white must be a logical vector.'; end

fdr = 0.05;
M = length(white);
K = sum(white);
GS = length(gsets.Name);

pval = nan(GS, 1);

for k = 1:GS
	gset = gsets.genes(k);
	if isempty(gset), continue, end
		
	valid(k) = true;
	
	white_in_gset = sum(white(gset));

	pval(k) = 1 - hygecdf(white_in_gset - 1, M, K, length(gset));
end

valid_gsets = find(~isnan(pval));
pval = pval(valid_gsets);

[~, order] = sort(pval, 'ascend');
valid_gsets = valid_gsets(order);
pval = pval(order);

num_significant = 0;
for k = length(pval):-1:1
	if pval(k) < k * fdr / length(pval)
		num_significant = k;
	end
end

num_significant


enriched = struct;
enriched.Name = gsets.Name(valid_gsets(1:num_significant));
for k = 1:num_significant
	enriched.Genes{k} = gsets.genes(valid_gsets(k));
end
enriched.Pval = pval(1:num_significant);

fprintf(1, 'Geneset\tNgenes\tP-value\n');
for k = 1:length(enriched.Name)
	fprintf(1, '%s\t%d\t%.8f\n', enriched.Name{k}, ...
		length(enriched.Genes{k}), enriched.Pval(k));
end





