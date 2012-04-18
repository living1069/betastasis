
function found_genes = vcfa_top_genes(vcfa_file, num_genes)

if nargin < 2
	num_genes = 100;
end

[data, headers] = readtable(vcfa_file);
genes = data{rx(headers, 'NEARBY_GENES')};

found_genes = {};

for k = 1:length(genes)
	cols = textscan(genes{k}, '%s', 'Delimiter', ';');
	for c = 1:length(cols)
		found_genes{end+1} = cols{c}{1};
	end
	if length(found_genes) > num_genes * 10, break, end
end

found_genes = unique_preserve_order(found_genes);
found_genes = found_genes(1:num_genes)';



