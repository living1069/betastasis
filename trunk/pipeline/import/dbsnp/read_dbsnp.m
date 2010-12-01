function dbsnp = read_dbsnp()

global organism;

chr_snps = repmat({struct}, length(organism.Chromosomes.Name), 1);

files = dir();
for k = 1:length(files)
	tokens = regexpi(files(k).name, '^ds_flat_ch(.+)\.flat$', 'tokens');
	if length(tokens) ~= 1, continue, end
	
	token = tokens{1}; chr = token{1};
	filepath = files(k).name;
	
	chr = chromosome_sym2num(chr);
	chr_snps{chr} = parse_dbsnp_chromosome(filepath, chr);
end

dbsnp = cat_structs(chr_snps{:});







function snp = parse_dbsnp_chromosome(filepath, chr)

snp_prealloc = 100000;

snp = struct;
snp.ID = cell(snp_prealloc, 1);
snp.Validated = false(snp_prealloc, 1);
snp.Chromosome = chr * ones(snp_prealloc, 1);
snp.Offset = zeros(snp_prealloc, 1);
snp.Alleles = cell(snp_prealloc, 1);

N = 0;

fid = fopen(filepath);
while 1
	line = fgetl(fid);
	if line == -1, break, end
	
	tokens = regexpi(line, '^(rs\d+) \| human', 'tokens');
	if length(tokens) == 1
		token = tokens{1}; id = token{1};
		N = N + 1;
		snp.ID{N} = id;
		continue;
	end
	
	tokens = regexpi(line, '^SNP \| alleles=''(.+)'' \|', 'tokens');
	if length(tokens) == 1
		token = tokens{1}; alleles = token{1};
		alleles = strrep(alleles, '/', '');
		snp.Alleles{N} = alleles;
		continue;
	end
	
	if regexpi(line, '^VAL \| validated=YES')
		snp.Validated(N) = true;
		continue;
	end
	
	tokens = regexpi(line, ...
		'^CTG \| .* \| chr-pos=(\d+) \| .* \| orient=(.)', 'tokens');
	if length(tokens) == 1
		token = tokens{1};
		offset = str2double(token{1}); strand = token{2};
		snp.Offset(N) = offset;
		if strcmp(strand, '-')
			snp.Alleles{N} = seqcomplement(snp.Alleles{N});
		elseif ~strcmp(strand, '+')
			error 'Unsupported strand.';
		end
		continue;
	end
end

fclose(fid);

[~, order] = sort(snp.Offset(1:N));

snp.ID = snp.ID(order);
snp.Validated = snp.Validated(order);
snp.Chromosome = snp.Chromosome(order);
snp.Offset = snp.Offset(order);
snp.Alleles = snp.Alleles(order);

