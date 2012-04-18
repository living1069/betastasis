function meth = import_tcga_level3_methylation(files)

if nargin < 1
	files = find_files('.*_methylation_analysis.txt');
end

meth.meta.tcga_barcode = {};
meth.meta.sample_id = {};
meth.beta = [];
meth.rows.gene_symbol = {};
meth.rows.chromosome = [];
meth.rows.position = [];

progress = Progress;

for f = 1:length(files)
	progress.update(f / length(files));
	
	[data, headers] = readtable(files{f});
	
	% There seems to be a weird artifact in TCGA data whereby some samples
	% have genes marked as DISCONTUNUED.
	valid = ~strcmp(data{3}, 'DISCONTUNUED') & ~strcmp(data{3}, '');
	for k = 1:length(data), data{k} = data{k}(valid); end

	samples = data{1};
	beta = str2double(data{2});
	gene = data{3};
	chromosome = chromosome_sym2num(data{4});
	position = str2double(data{5});

	[~, sample_starts] = unique(samples, 'first');
	S = length(sample_starts);
	sample_starts = [sample_starts; length(samples)+1];
	
	for s = 1:length(sample_starts)-1
		meth.meta.tcga_barcode{end+1} = samples{s};
		meth.meta.sample_id{end+1} = samples{s}(1:15);
		
		lines = sample_starts(s):sample_starts(s+1)-1;
		if isempty(meth.rows.gene_symbol)
			meth.rows.gene_symbol = gene(lines);
			meth.rows.chromosome = chromosome(lines);
			meth.rows.position = position(lines);
		else
			%meth.rows
			%whos gene
			%meth.rows.gene_symbol(end-29:end)
			%gene(end-29:end)
			
			% Check that gene names match between samples.
			if any(~strcmp(meth.rows.gene_symbol, gene(lines)))
				error 'Gene names do not match between samples.';
			end
		end
		
		meth.beta(:, length(meth.meta.sample_id)) = beta(lines);
	end
end

meth.meta.type = 'Methylation microarray beta-values';



