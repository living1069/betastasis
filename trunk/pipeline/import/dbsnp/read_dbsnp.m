function [] = read_dbsnp()

files = dir();
for k = 1:length(files)
	tokens = regexpi(files(k).name, '^ds_flat_ch(.+)\.flat$', 'tokens');
	if length(tokens) ~= 1, continue, end
	
	token = tokens{1}; chr = token{1};
	filepath = files(k).name;
	
	fprintf(1, 'Reading chromosomes from file %s...\n', filepath);
	chr = chromosome_sym2num(chr);
	parse_dbsnp_chromosome(filepath, chr);
end








function [] = parse_dbsnp_chromosome(filepath, chr)

global organism;

snp_prealloc = 10^7;

snps = struct;
snps.ID = zeros(snp_prealloc, 1);
snps.Validated = false(snp_prealloc, 1);
snps.Offset = zeros(snp_prealloc, 1);
snps.Alleles = cell(snp_prealloc, 1);

N = 0;

fid = fopen(filepath);

while 1
	line = fgetl(fid);
	if line == -1, break, end
		
	if length(line) < 4 || strcmp(line(1:2), 'ss'), continue, end
	
	if strcmp(line(1:2), 'rs')
		tokens = regexpi(line, '^rs(\d+) \|', 'tokens');
		if length(tokens) == 1
			token = tokens{1}; id = str2double(token{1});
			N = N + 1;
			snps.ID(N) = id;
		end
		continue;
	end
	
	code = line(1:3);
	
	if strcmp(code, 'SNP')
		tokens = regexpi(line, '^SNP \| alleles=''(.+)'' \|', 'tokens');
		if length(tokens) == 1
			token = tokens{1}; alleles = token{1};
			snps.Alleles{N} = alleles;
		end
		continue;
	end
	
	if strcmp(code, 'CTG')
		tokens = regexpi(line, '^CTG \| assembly=GRCh37 .* \| chr-pos=(\d+)',...
			'tokens');
		if length(tokens) == 1
			token = tokens{1};
			snps.Offset(N) = str2double(token{1});
		end
		continue;
	end
	
	if regexpi(line, '^VAL \| validated=YES')
		snps.Validated(N) = true;
		continue;
	end
end

fclose(fid);

[~, order] = sort(snps.Offset(1:N));

snps.ID = snps.ID(order);
snps.Validated = snps.Validated(order);
snps.Offset = snps.Offset(order);
snps.Alleles = snps.Alleles(order);

save(sprintf('chr%s', organism.Chromosomes.Name{chr}), 'snps');

