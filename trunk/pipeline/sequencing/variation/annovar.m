
% Author: Matti Annala <matti.annala@tut.fi>

function [] = annovar(variant_file)

global organism;
chromosomes = organism.Chromosomes;

if ~rx(variant_file, '\.vcf')
	error 'Only VCF files are currently supported.';
end

tmp = temporary('annovar');
out_prefix = [tmp 'variants']
	
[~, headers] = readtable(variant_file, 'Comment', '^##', 'NumLines', 0);
samples = headers((find(strcmp(headers, 'FORMAT'))+1):end);
samples = regexprep(samples, '.*/', '');

% Command used to convert the VCF file into an ANNOVAR compatible format.
% Note that the 'sed' command is used to fix a bug in convert2annovar.pl
convert_cmd = sprintf(['convert2annovar.pl --format vcf4 --includeinfo ' ...
	'--allallele <(gunzip -c %s) | sed "s/\t\t/\t/"'], variant_file);

humandb_dir = '~/tools/annovar*/humandb';
	
% First we calculate genomic annotations.
unix(sprintf([ ...
	'annotate_variation.pl --buildver hg19 --geneanno --outfile %s ' ...
	'<(%s) %s'], out_prefix, convert_cmd, humandb_dir));

% Determine which variants were also found by the 1000 Genomes project.
unix(sprintf([ ...
	'annotate_variation.pl --buildver hg19 --filter --outfile %s ' ...
	'--dbtype 1000g2011may_all <(%s) %s'], ...
	out_prefix, convert_cmd, humandb_dir));

% Determine which variants were also found by the dbSNP project.
unix(sprintf([ ...
	'annotate_variation.pl --buildver hg19 --filter --outfile %s ' ...
	'--dbtype snp132 <(%s) %s'], out_prefix, convert_cmd, humandb_dir));
	
% Determine which variants are reported in COSMIC.
unix(sprintf([ ...
	'annotate_variation.pl --buildver hg19 --regionanno --outfile %s ' ...
	'--dbtype bed -bedfile cosmic_v56.bed <(%s) %s'], ...
	[out_prefix '.cosmic'], convert_cmd, humandb_dir));
	

	

% Construct mappings from genomic loci to annotated features.
data = readtable([out_prefix '.exonic_variant_function'], ...
	'Header', false, 'Ignore', [1 9:1000]);
exonic_func_map = containers.Map( ...
	strcat(data{3},':',data{4},'-',data{5},':',data{6},'>',data{7}), ...
	strcat(data{1}, {sprintf('\t')}, data{2}));

data = readtable([out_prefix '.hg19_snp132_dropped'], ...
	'Header', false, 'Ignore', [1 8:1000]);
dbsnp_map = containers.Map( ...
	strcat(data{2},':',data{3},'-',data{4},':',data{5},'>',data{6}), ...
	data{1});

data = readtable([out_prefix '.hg19_ALL.sites.2011_05_dropped'], ...
	'Header', false, 'Ignore', [1 8:1000]);
kgenomes_map = containers.Map( ...
	strcat(data{2},':',data{3},'-',data{4},':',data{5},'>',data{6}), ...
	data{1});
	
data = readtable([out_prefix '.cosmic.hg19_bed'], 'Header', false, ...
	'Ignore', [1 2 6:1000]);
cosmic_map = containers.Map(strcat(data{1},':',data{2},'-',data{3}), ...
	repmat({'YES'}, 1, length(data{1})));

	
	
	
	
	
	
	
% Finally, we combine all of the lists into one gigantic pile of variant info.
data = readtable([out_prefix '.variant_function']);
fid = fopen('variants.vcfa', 'W');
fprintf(fid, ['CHROMOSOME\tSTART\tEND\tREFERENCE\tOBSERVED\t' ...
	'FUNCTION\tNEARBY_GENES\tSYNONYMOUS\tPROTEIN_EFFECT\tDETAILS\t' ...
	'DBSNP\t1000GENOMES\tCOSMIC']);
fprintf(fid, '\t%s', samples{:});
fprintf(fid, '\n');
	
keys = strcat(data{3},':',data{4},'-',data{5},':',data{6},'>',data{7});

valid = exonic_func_map.isKey(keys);
exonic_funcs = repmat({sprintf('-\t-')}, length(data{1}), 1);
exonic_funcs(valid) = exonic_func_map.values(keys(valid));

valid = dbsnp_map.isKey(keys);
dbsnps = repmat({'-'}, length(data{1}), 1);
dbsnps(valid) = dbsnp_map.values(keys(valid));
fprintf('%d / %d (%.1f%%) variants found in dbSNP v132.\n', ...
	sum(valid), length(valid), sum(valid) / length(valid) * 100);

valid = kgenomes_map.isKey(keys);
kgenomes = repmat({'-'}, length(data{1}), 1);
kgenomes(valid) = kgenomes_map.values(keys(valid));
fprintf('%d / %d (%.1f%%) variants found in 1000 Genomes.\n', ...
	sum(valid), length(valid), sum(valid) / length(valid) * 100);
	
% Region based annotations: nucleotides don't matter, only position does.
keys = strcat(data{3},':',data{4},'-',data{5});

valid = cosmic_map.isKey(keys);
cosmic = repmat({'-'}, length(data{1}), 1);
cosmic(valid) = cosmic_map.values(keys(valid));
fprintf('%d / %d (%.1f%%) variants found in COSMIC.\n', ...
	sum(valid), length(valid), sum(valid) / length(valid) * 100);

	
% Now we calculate the reference sequence surrounding the variant site.
%chr = chromosome_sym2num(data{3});
%pos = [str2double(data{4}), str2double(data{5})];
%flank_seq = cell(length(chr), 1);
%for k = 1:length(flank_seq)
%	flank_seq{k} = lower(chromosomes.Sequence{chr(k)}(pos(k,1)-30:pos(k,2)+30));
%	if strcmp(data{6}{k}, '-')
%		flank_seq{k} = [flank_seq{k}(2:31) '-' flank_seq{k}(32:end)];
%	else
%		flank_seq{k}(31:31+pos(k,2)-pos(k,1)) = ...
%			upper(flank_seq{k}(31:31+pos(k,2)-pos(k,1)));
%	end
%end

% Double check that we are picking the sample data from the right column.
for k = 1:length(data)
	if rx(data{k}{1}, '[01]/[01]:.*:')
		sample_col = k; break;
	end
end

for k = 1:length(data{1})
	fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
		data{3}{k}, data{4}{k}, data{5}{k}, data{6}{k}, data{7}{k}, ...
		data{1}{k}, data{2}{k}, exonic_funcs{k}, data{15}{k}, dbsnps{k}, ...
		kgenomes{k}, cosmic{k});

	for c = [sample_col:length(data)]
		fprintf(fid, '\t%s', data{c}{k});
	end
	fprintf(fid, '\n');
end

fclose(fid);



