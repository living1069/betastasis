function merged = merge_bed_features(bed_file, method)

data = readtable(bed_file, 'Numeric', 2:3);

[~, order] = sort(data{4});

for k = 1:length(data), data{k} = data{k}(order); end

% Throw away nonsensical chromosomes.
chr = chromosome_sym2num(data{1});
valid = ~isnan(chr);
for k = 1:length(data), data{k} = data{k}(valid); end
chr = chr(valid);
pos = [data{2} data{3}];
sorted_genes = data{4};

run_ends = [ find(~strcmp(sorted_genes(1:end-1), sorted_genes(2:end))); ...
	length(sorted_genes) ];
run_starts = [1; run_ends(1:end-1)+1];

if rx(method, 'union')
else
	error 'Unsupported method requested for merging BED features.';
end

E = 0;

merged = struct;
merged.chromosome = nan(length(chr), 1);
merged.position = nan(length(chr), 2);
merged.feature = cell(length(chr), 1);

for r = [run_starts'; run_ends']
	gene_exons = nan(0, 3);   % [chr, start, end] matrix
	
	for k = r(1):r(2)
		% For each exon, we check if it overlaps with an existing exon
		overlapping = find(gene_exons(:, 1) == chr(k) & ...
			gene_exons(:, 2) < pos(k, 2) & gene_exons(:, 3) > pos(k, 1));
		if isempty(overlapping)
			gene_exons(end+1, :) = [chr(k) pos(k, :)];
		else
			gene_exons(end+1, :) = [chr(k) ...
				min([gene_exons(overlapping, 2); pos(k, 1)]) ...
				max([gene_exons(overlapping, 3); pos(k, 2)])];
			gene_exons(overlapping, :) = [];
		end
	end
	
	En = size(gene_exons, 1);
	merged.chromosome(E+1:E+En) = gene_exons(:, 1);
	merged.position(E+1:E+En, :) = gene_exons(:, 2:3);
	merged.feature(E+1:E+En) = repmat({sorted_genes{r(1)}}, En, 1);
	E = E + En;
end

merged = filter_rows(merged, ~isnan(merged.chromosome));




% SOME BASH CODE FOR CONVERTING ENSEMBL GTF TO A SUITABLE BED FILE
%grep '\sexon\s' Homo_sapiens.GRCh37.67.gtf | cut -f 1,4,5,9 | sed 's/gene_id "\([^\s]*\)".*/\1/'

% UNRELATED
%paste <(gunzip -c variants.vcfa.gz | tail -n +2 | cut -f 1-8)
%<(gunzip -c variants.vcfa.gz | tail -n +2 | cut -f 9-)
