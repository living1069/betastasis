
function variants = read_crappy_vcf(vcf_file)

[data, headers] = readtable(vcf_file, 'Comment', '^##', ...
	'Ignore', '^(ID|QUAL|FILTER|INFO|FORMAT)$', 'Numeric', '^POSITION$');

for k = 1:length(headers)
	if ~iscellstr(data{k}), continue, end
	if any(rx(data{k}(1:100), '[01]/[01](:| \()'))
		first_sample_col = k;
		break;
	end
end

S = length(headers) - first_sample_col + 1;
V = length(data{1});

variants = struct;
variants.meta.sample_id = headers(first_sample_col:end)';
variants.rows.chromosome = chromosome_sym2num(data{rx(headers, '^CHROM')});
variants.rows.position = data{rx(headers, '^(POSITION|POS)$')};
variants.rows.ref_allele = data{rx(headers, '^(REFERENCE|REF)$')};
variants.rows.alt_allele = data{rx(headers, '^(ALTERNATE|ALT)$')};

if any(rx(headers, '^FUNCTION'))
	variants.rows.function = data{rx(headers, '^FUNCTION')};
end
if any(rx(headers, 'NEARBY.*GENES'))
	variants.rows.nearby_genes = data{rx(headers, 'NEARBY.*GENES')};
end
if any(rx(headers, 'EXONIC_FUNCTION'))
	variants.rows.synonymous = data{rx(headers, 'EXONIC_FUNCTION')};
end
if any(rx(headers, 'AA_CHANGE'))
	variants.rows.protein_effect = data{rx(headers, 'AA_CHANGE')};
end
if any(rx(headers, 'dbSNP'))
	variants.rows.dbsnp = data{rx(headers, 'dbSNP')};
end
if any(rx(headers, '1000G'))
	variants.rows.kgenomes = data{rx(headers, '1000G')};
end
if any(rx(headers, 'ESP6500'))
	variants.rows.esp6500 = data{rx(headers, 'ESP6500')};
end
if any(rx(headers, 'COSMIC'))
	variants.rows.cosmic = data{rx(headers, 'COSMIC')};
end

variants.genotype = nan(V, S);
%variants.genotype_quality = nan(V, S);



for s = 1:S
	info = data{first_sample_col + s - 1};
	%if any(info{1} == ':')
		str = sprintf('%s\n', info{:});
		sdata = textscan(str, '%s', 'Delimiter', '\n', 'ReturnOnError', false);
		gtypes = char(sdata{1});
		valid = any(gtypes == '/', 2);
		variants.genotype(valid, s) = sum(gtypes(valid, 1:3) == '1', 2);
		
		%if s == 1, str, sdata{1}, gtypes, end

		%variants.genotype_quality(:, s) = max(sdata{2:3});
	%else
	%	valid = ~cellfun(@isempty, info);
	%	tokens = rx_capture(info(valid), '^(./.) \((\d+)\)$');
	%	variants.genotype(valid, s) = sum(char(tokens(:, 1)) == '1', 2);
	%	variants.genotype_quality(valid, s) = str2double(tokens(:, 2));
	%end
end

