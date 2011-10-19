
% Author: Matti Annala <matti.annala@tut.fi>

function variants = import_vcf(vcf_file)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;
exons = organism.Exons;

variant_phred_threshold = 30;

fid = fopen(vcf_file);
while 1
	line = fgetl(fid);
	if ~ischar(line), error 'Abrupt file termination.', end
	if ~strcmp(line(1:2), '##'), break, end
end

cols = textscan(line(2:end), '%s', 'Delimiter', '\t');
headers = cols{1};

samples = headers(10:end);
samples = regexprep(samples, '^.*\/(.+).bam$', '$1');
S = length(samples);

format_str = '';
for k = 1:length(headers), format_str = [format_str '%s']; end

	
	
data = textscan(fid, format_str, 'Delimiter', '\t', 'BufSize', 16384, ...
	'ReturnOnError', false);
fclose(fid);


% Filter out variants with a non-sufficient quality.
valid = str2double(data{6}) > variant_phred_threshold;
for k = 1:length(data), data{k} = data{k}(valid); end
V = sum(valid);

chr = chromosome_sym2num(data{1});
pos = str2double(data{2});


variants = struct;
variants.Samples = samples';
variants.Chromosome = chr;
variants.Position = pos;
variants.Alleles = cell(V, 1);
variants.OverlapGenes = cell(V, 1);
variants.ReferenceSeq = cell(V, 1);
variants.Consequences = cell(V, 1);
variants.Genotype = nan(V, S);
variants.TotalReads = nan(V, S);







sample_data = data(10:end);
for s = 1:S
	
	tokens = regexp(sample_data{s}, '(.+?):.+?:(\d+):\d+', 'tokens');
	
	gtype = cell(length(tokens), 1);
	read_depth = nan(length(tokens), 1);
	for k = 1:length(tokens)
		token = tokens{k}; token = token{1};
		gtype{k} = token{1}; read_depth(k) = str2double(token{2});
	end
		
	variants.Genotype(strcmp(gtype, '0/0'), s) = 0;
	variants.Genotype(strcmp(gtype, '1/0') | strcmp(gtype, '0/1'), s) = 1;
	variants.Genotype(strcmp(gtype, '1/1'), s) = 2;
	
	%sdata = textscan(sample_data{s}, '%s %*d %d %*d', 'Delimiter', ':');
		
	variants.TotalReads(:, s) = read_depth;
end



% Calculate which genes overlap with the SNV, and what the consequences are.
genes_in_chr = {};
for c = 1:length(chromosomes.Name)
	genes_in_chr{c} = find(genes.Chromosome == c);
end

inside_gene = cell(V, 1);
for k = 1:V
	variants.Alleles{k} = sprintf('%s>%s', data{4}{k}, data{5}{k});
	
	inside_gene{k} = find(genes.Position(genes_in_chr{chr(k)}, 1) <= pos(k) &...
		genes.Position(genes_in_chr{chr(k)}, 2) >= pos(k))';
	variants.OverlapGenes{k} = genes_in_chr{chr(k)}(inside_gene{k});
	
	seq = lower(chromosomes.Sequence{chr(k)}(pos(k)-20:pos(k)+20));
	seq(21) = upper(seq(21));
	variants.ReferenceSeq{k} = seq;
	
	variants.Consequences{k} = snv_consequence(sprintf('chr%s:%d:%s', ...
		chromosomes.Name{chr(k)}, pos(k), variants.Alleles{k}));
end




