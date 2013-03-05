
% Author: Matti Annala <matti.annala@tut.fi>

function [] = annovar(variant_file)

if ~rx(variant_file, '\.vcf')
	error 'Only VCF files are currently supported.';
end

tmp = temporary('annovar');
out_prefix = [tmp 'variants']
	
[~, headers] = readtable(variant_file, 'Comment', '^##', 'NumLines', 0);

% Sanity check: does the file actually look like it's in VCF format?
if any(~strcmp({'#CHROM', 'POS', 'ID', 'REF', 'ALT'}', headers(1:5)))
	error 'File does not look like a VCF file.';
end

extra_col = 6;
headers(extra_col:end) = regexprep(headers(extra_col:end), ...
	'(.*/)?(.*).bam', '$2');

if rx(variant_file, '\.gz$')
	variant_file = sprintf('<(gunzip -c %s)', variant_file);
end

% Command used to convert the VCF file into an ANNOVAR compatible format.
% Note that the 'sed' command is used to fix a bug in convert2annovar.pl
convert_cmd = sprintf(['convert2annovar.pl --format vcf4 --includeinfo ' ...
	'--allallele %s | sed "s/\t\t/\t/"'], variant_file);

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
	
% Determine which variants were also found by the ESP6500 project.
unix(sprintf([ ...
	'annotate_variation.pl --buildver hg19 --filter --outfile %s ' ...
	'--dbtype esp6500si_all <(%s) %s'], out_prefix, convert_cmd, humandb_dir));
	
% Determine which variants are reported in COSMIC.
unix(sprintf([ ...
	'annotate_variation.pl --buildver hg19 --regionanno --outfile %s ' ...
	'--dbtype bed -bedfile cosmic_v56.bed <(%s) %s'], ...
	[out_prefix '.cosmic'], convert_cmd, humandb_dir));

% Construct mappings from genomic loci to annotated features.
fprintf('Reading exonic variant annotations...\n');
[out_prefix '.exonic_variant_function']
data = readtable([out_prefix '.exonic_variant_function'], ...
	'Header', false, 'Ignore', [1 9:1000]);
exonic_func_keys = strcat(data{3},':',data{4},':',data{6},'>',data{7})
exonic_func_values = strcat(data{1}, {sprintf('\t')}, data{2});

fprintf('Reading dbSNP annotations...\n');
data = readtable([out_prefix '.hg19_snp132_dropped'], ...
	'Header', false, 'Ignore', [1 8:1000]);
dbsnp_keys = strcat(data{2},':',data{3},':',data{5},'>',data{6});
dbsnp_values = data{1};

fprintf('Reading 1000 Genomes annotations...\n');
data = readtable([out_prefix '.hg19_ALL.sites.2011_05_dropped'], ...
	'Header', false, 'Ignore', [1 8:1000]);
kgenomes_keys = strcat(data{2},':',data{3},':',data{5},'>',data{6});kgenomes_values = data{1};

fprintf('Reading ESP6500 annotations...\n');
data = readtable([out_prefix '.hg19_esp6500si_all_dropped'], ...
	'Header', false, 'Ignore', [1 8:1000]);
esp6500_keys = strcat(data{2}, ':', data{3}, ':', data{5}, '>', data{6});
esp6500_values = data{1};
	
fprintf('Reading COSMIC annotations...\n');
data = readtable([out_prefix '.cosmic.hg19_bed'], 'Header', false, ...
	'Ignore', [1 2 6:1000]);
cosmic_keys = strcat(data{1},':',data{2});
cosmic_values = repmat({'YES'}, 1, length(data{1}));

	
	
	
	
	
	
	
% Finally, we combine all of the lists into one gigantic pile of variant info.
fprintf('Reading main ANNOVAR variant information...\n');
data = readtable([out_prefix '.variant_function'], ...
	'Header', false, 'Ignore', [5 8:11]);
	
variant_func = data(1:2);
data = data(3:end);
	
keys = strcat(data{1},':',data{2},':',data{3},'>',data{4});

[valid, pos] = ismember(keys, exonic_func_keys);
exonic_funcs = repmat({sprintf('-\t-')}, length(keys), 1);
exonic_funcs(valid) = exonic_func_values(pos(valid));

[valid, pos] = ismember(keys, dbsnp_keys);
dbsnps = repmat({'-'}, length(data{1}), 1);
dbsnps(valid) = dbsnp_values(pos(valid));
fprintf('%d / %d (%.1f%%) variants found in dbSNP v132.\n', ...
	sum(valid), length(valid), sum(valid) / length(valid) * 100);

[valid, pos] = ismember(keys, kgenomes_keys);
kgenomes = repmat({'-'}, length(data{1}), 1);
kgenomes(valid) = kgenomes_values(pos(valid));
fprintf('%d / %d (%.1f%%) variants found in 1000 Genomes.\n', ...
	sum(valid), length(valid), sum(valid) / length(valid) * 100);
	
[valid, pos] = ismember(keys, esp6500_keys);
esp6500 = repmat({'-'}, length(data{1}), 1);
esp6500(valid) = esp6500_values(pos(valid));
fprintf('%d / %d (%.1f%%) variants found in ESP6500.\n', ...
	sum(valid), length(valid), sum(valid) / length(valid) * 100);
	
% Region based annotations: nucleotides don't matter, only position does.
keys = strcat(data{1},':',data{2});

[valid, pos] = ismember(keys, cosmic_keys);
cosmic = repmat({'-'}, length(data{1}), 1);
cosmic(valid) = cosmic_values(pos(valid));
fprintf('%d / %d (%.1f%%) variants found in COSMIC.\n', ...
	sum(valid), length(valid), sum(valid) / length(valid) * 100);

	
	
fid = fopen('variants.vcfa', 'W');
fprintf(fid, ['CHROM\tPOS\tID\tREF\tALT\t' ...
	'FUNCTION\tNEARBY_GENES\tSYNONYMOUS\tPROTEIN_EFFECT\t' ...
	'DBSNP\t1000GENOMES\tESP6500\tCOSMIC']);
fprintf(fid, '\t%s', headers{extra_col:end});
fprintf(fid, '\n');
	
for k = 1:length(data{1})
	fprintf(fid, '%s\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
		data{1}{k}, data{2}{k}, data{3}{k}, data{4}{k}, ...
		variant_func{1}{k}, variant_func{2}{k}, exonic_funcs{k}, dbsnps{k}, ...
		kgenomes{k}, esp6500{k}, cosmic{k});

	% For some reason ANNOVAR removes the "ID" column so that's why -1
	for c = [extra_col:length(data)]
		fprintf(fid, '\t%s', data{c}{k});
	end
	fprintf(fid, '\n');
end

fclose(fid);



