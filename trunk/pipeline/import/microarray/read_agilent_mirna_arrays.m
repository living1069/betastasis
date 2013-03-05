function mirna_expr = read_agilent_mirna_arrays(sample_files)

S = length(sample_files);

mirna_expr = struct;
mirna_expr.meta.sample_id = sample_files(:)';
mirna_expr.mean = nan(0, S);

for s = 1:S
	[data, headers] = readtable(sample_files{s}, 'HeaderRegex', '^FEATURES', ...
		'IncludeColumns', {'SystematicName', 'ProcessedSignal'}, ...
		'Numeric', 'ProcessedSignal');

	if s == 1
		mirna_expr.rows.gene = data{rx(headers, 'SystematicName')};
		G = length(mirna_expr.rows.gene);
		mirna_expr.mean(G, S) = 0;   % Preallocate space
	else
		if ~strcmp(data{rx(headers, 'SystematicName')}, mirna_expr.rows.gene)
			error('MicroRNA annotations do not match between samples.');
		end
	end
	
	mirna_expr.mean(:, s) = data{rx(headers, 'ProcessedSignal')};
end

% Discard dark corner probes etc.
mirna_expr = filter_rows(mirna_expr, ...
	~rx(mirna_expr.rows.gene, 'Corner|Control'));


