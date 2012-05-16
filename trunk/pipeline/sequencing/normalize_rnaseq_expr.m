
function norm = normalize_rnaseq_expr(expr, method)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;

if nargin < 2, method = 'rpkm'; end

S = size(expr.mean, 2);

norm = expr;

if rx(method, 'PKM')
	norm.meta.expr_normalization = repmat({ 'RPKM' }, 1, S);
	gidx = gene_idx(expr.rows.gene_symbol);
	valid = ~isnan(gidx);
	if any(~valid)
		norm = filter_rows(norm, valid);
		gidx = gidx(valid);
		fprintf('Discarded %d genes with an unknown name.\n', sum(~valid));
	end
	
	gene_tx_len = nan(length(gidx), 1);
	for k = 1:length(gidx)
		gene_tx_len(k) = mean(cellfun(@length, ...
			transcripts.Sequence(transcripts.Gene == gidx(k))));
	end
	
	for s = 1:S
		total_aligned = sum(expr.mean(:, s));
		norm.mean(:, s) = expr.mean(:, s) ./ (gene_tx_len / 1000) / ...
			(total_aligned / 1e6);
	end
else
	error 'Unsupported normalization method requested.';
end


