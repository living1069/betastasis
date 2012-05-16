
function variants = read_vcf(vcf_file)

[data, headers] = readtable(vcf_file, 'Comment', '^##', ...
	'Ignore', '^(ID|QUAL|FILTER|FORMAT)$', 'Numeric', '^(START|POS)$');

for k = 1:length(headers)
	if ~iscellstr(data{k}), continue, end
	if rx(data{k}{1}, '[01]/[01]:')
		first_sample_col = k;
		break;
	end
end

S = length(headers) - first_sample_col + 1;
V = length(data{1});

variants = struct;
variants.meta.sample_id = headers(first_sample_col:end)';
variants.rows.chromosome = chromosome_sym2num(data{rx(headers, '^#?CHROM')});
variants.rows.position = data{rx(headers, '^(START|POS)$')};
variants.rows.ref_allele = data{rx(headers, '^(REFERENCE|REF)$')};
variants.rows.alt_allele = data{rx(headers, '^(OBSERVED|ALT)$')};

if any(rx(headers, 'FLANK.*SEQ'))
	variants.rows.flanking_seq = data{rx(headers, 'FLANK.*SEQ')};
end
if any(rx(headers, 'FUNCTION'))
	variants.rows.function = data{rx(headers, 'FUNCTION')};
end
if any(rx(headers, 'NEARBY.*GENES'))
	variants.rows.nearby_genes = data{rx(headers, 'NEARBY.*GENES')};
end
if any(rx(headers, 'SYNONYMOUS'))
	variants.rows.synonymous = data{rx(headers, 'SYNONYMOUS')};
end
if any(rx(headers, 'PROTEIN.*EFFECT'))
	variants.rows.protein_effect = data{rx(headers, 'PROTEIN.*EFFECT')};
end
%variants.rows.kgenomes = ~strcmpi(data{rx(headers, '1000GENOMES')}, '-');
%variants.rows.cosmic = ~strcmpi(data{rx(headers, 'COSMIC')}, '-');

variants.genotype = nan(V, S);
variants.genotype_quality = nan(V, S);

for s = 1:S
	info = data{first_sample_col + s - 1};
	str = sprintf('%s\n', info{:});
	sdata = textscan(str, '%s %d %d %d %*[^\n]', 'Delimiter', ':,', ...
		'ReturnOnError', false);
	
	variants.genotype(:, s) = sum(char(sdata{1}) == '1', 2);
	variants.genotype_quality(:, s) = max(sdata{2}, sdata{3});
end

