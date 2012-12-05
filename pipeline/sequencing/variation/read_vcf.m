
function variants = read_vcf(vcf_file)

[data, headers] = readtable(vcf_file, 'Comment', '^##', ...
	'Ignore', '^(ID|QUAL|FILTER|INFO|FORMAT)$', 'Numeric', '^(START|POS)$');

if any(rx(headers, 'HOM/HET/TOT'))    % Genotype style = 0 / 1 (255)
	first_sample_col = find(rx(headers, 'HOM/HET/TOT'), 1, 'last') + 1
else                             % Genotype style = 0/1:...
	for k = 1:length(headers)
		if ~iscellstr(data{k}), continue, end
		if any(rx(data{k}(1:10), '[01]/[01][: ]'))
			first_sample_col = k;
			break;
		end
	end
end

S = length(headers) - first_sample_col + 1;
V = length(data{1});

variants = struct;
variants.meta.sample_id = headers(first_sample_col:end)';
variants.rows.chromosome = chromosome_sym2num(data{rx(headers, '^#*CHR')});
variants.rows.position = data{rx(headers, '^(START|POS|POSITION)$')};
variants.rows.ref_allele = data{rx(headers, '^(REFERENCE|REF)$')};
variants.rows.alt_allele = data{rx(headers, '^(OBSERVED|OBS|ALT)$')};

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
if any(rx(headers, 'DBSNP'))
	variants.rows.dbsnp = data{rx(headers, 'DBSNP')};
end
if any(rx(headers, '1000GENOMES'))
	variants.rows.kgenomes = str2double(data{rx(headers, '1000GENOMES')});
	if all(isnan(variants.rows.kgenomes))   % FIXME
		variants.rows.kgenomes = ...
			strcmp(data{rx(headers, '1000GENOMES')}, 'YES');
	end
end
if any(rx(headers, 'COSMIC'))
	variants.rows.cosmic = ~strcmpi(data{rx(headers, 'COSMIC')}, '-');
end

variants.genotype = nan(V, S);
variants.genotype_quality = nan(V, S);

for s = 1:S
	info = data{first_sample_col + s - 1};
	if any(info{1} == ':')
		str = sprintf('%s\n', info{:});
		sdata = textscan(str, '%s %d %d %d %*[^\n]', 'Delimiter', ':,', ...
			'ReturnOnError', false);
		variants.genotype(:, s) = sum(char(sdata{1}) == '1', 2);
		variants.genotype_quality(:, s) = max(sdata{2:3});
	else
		valid = ~cellfun(@isempty, info);
		tokens = rx_capture(info(valid), '^(./.) \((\d+)\)$');
		variants.genotype(valid, s) = sum(char(tokens(:, 1)) == '1', 2);
		variants.genotype_quality(valid, s) = str2double(tokens(:, 2));
	end
end

