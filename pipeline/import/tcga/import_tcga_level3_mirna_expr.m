function expr = import_tcga_level3_mirna_expr(files)

if nargin < 1
	files = find_files('.*_mirna_expression_analysis.txt');
end

expr.meta.type = 'miRNA expression';
expr.meta.tcga_barcode = {};
expr.meta.sample_id = {};
expr.mean = [];
expr.rows.mirna_symbol = {};
expr.scale = {};

for f = 1:length(files)
	fid = fopen(files{f});
	data = textscan(fid, '%s %s %f', 'Headerlines', 1, ...
		'Delimiter', '\t', 'ReturnOnError', 0, 'TreatAsEmpty', 'NA');
	fclose(fid);
	
	samples = data{1};
	mirna = data{2};
	log_expr = data{3};

	[~, sample_starts] = unique(samples, 'first');
	S = length(sample_starts);
	sample_starts = [sample_starts; length(samples)+1];
	
	for s = 1:length(sample_starts)-1
		expr.meta.tcga_barcode{end+1} = samples{s};
		expr.meta.sample_id{end+1} = samples{s}(1:15);
		expr.scale{s} = 'Log-2';
		
		lines = sample_starts(s):sample_starts(s+1)-1;
		if isempty(expr.rows.mirna_symbol)
			expr.rows.mirna_symbol = mirna(lines);
		else
			% Check that gene names match between samples.
			if any(~strcmp(expr.rows.mirna_symbol, mirna(lines)))
				error 'MicroRNA names don''t match between samples.';
			end
		end
		
		expr.mean(:, length(expr.meta.sample_id)) = log_expr(lines);
	end
end

