function gene_expr = read_agilent_gene_arrays(sample_files)

S = length(sample_files);

gene_expr = struct;
gene_expr.rows = struct;
gene_expr.meta.sample_id = cell(1, S*2);
gene_expr.meta.sample_id(1:2:end) = strcat(sample_files(:), '(Cy5)');
gene_expr.meta.sample_id(2:2:end) = strcat(sample_files(:), '(Cy3)');
gene_expr.mean = nan(0, S*2);

progress = Progress;

for s = 1:S
	[data, headers] = readtable(sample_files{s}, 'HeaderRegex', '^FEATURES', ...
		'IncludeColumns', {'SystematicName', 'ProcessedSignal'}, ...
		'Numeric', 'ProcessedSignal');

	if s == 1
		gene_expr.rows.name = data{rx(headers, 'SystematicName')};
		G = length(gene_expr.rows.name);
		gene_expr.mean(G, S*2) = 0;   % Preallocate space
	else
		if ~strcmp(data{rx(headers, 'SystematicName')}, gene_expr.rows.name)
			error('Gene annotations do not match between samples.');
		end
	end
	
	gene_expr.mean(:, s*2-1) = data{rx(headers, 'rProcessedSignal')};
	gene_expr.mean(:, s*2  ) = data{rx(headers, 'gProcessedSignal')};
	
	progress.update(s / S);
end

% Discard dark corner probes etc.
gene_expr = filter_rows(gene_expr, ...
	~rx(gene_expr.rows.name, 'Corner|Control'));
gene_expr = filter_rows(gene_expr, ...
	rx(gene_expr.rows.name, '^(NM_|XM_|NR_|XR_|ENST)'));


