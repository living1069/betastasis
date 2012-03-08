function expr = import_tcga_level3_gene_expr(files)

if nargin < 1
	files = find_files('.*_gene_expression_analysis.txt');
end

expr.meta.sample_id = {};
expr.mean = [];
expr.rows.gene_symbol = {};

for f = 1:length(files)
	[data, headers] = readtable(files{f});
	samples = data{1};
	gene = data{2};
	log_expr = str2double(data{3});

	[~, sample_starts] = unique(samples, 'first');
	S = length(sample_starts);
	sample_starts = [sample_starts; length(samples)+1];
	
	for s = 1:length(sample_starts)-1
		expr.meta.sample_id{end+1} = samples{s};
		
		lines = sample_starts(s):sample_starts(s+1)-1;
		if isempty(expr.rows.gene_symbol)
			expr.rows.gene_symbol = gene(lines);
		else
			% Check that gene names match between samples.
			if any(~strcmp(expr.rows.gene_symbol, gene(lines)))
				error 'Gene names don''t match between samples.';
			end
		end
		
		expr.mean(:, length(expr.meta.sample_id)) = log_expr(lines);
	end
end

if sum(sum(expr.mean < 0)) / numel(expr.mean) > 0.1
	expr.meta.type = 'Differential gene expression';
else
	expr.meta.type = 'Gene expression';
end

expr.scale = repmat({'Natural'}, 1, size(expr.mean, 2));
expr.mean = 2.^expr.mean;

