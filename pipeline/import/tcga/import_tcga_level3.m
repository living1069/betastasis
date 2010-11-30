function [] = import_tcga_level3(filename, dataset)

if regexpi(filename, 'HumanMethylation27__methylation_analysis.txt')
	import_illumina_infinium(filename, dataset);
elseif regexpi(filename, 'OMA002_CPI__beta-value.txt')
	import_illumina_oma(filename, dataset, 2);
elseif regexpi(filename, 'OMA003_CPI__beta-value.txt')
	import_illumina_oma(filename, dataset, 3);
elseif regexpi(filename, '.*\.maf')
	import_tcga_mutations_maf(filename, dataset);
else
	error 'Could not recognize the type of the level 3 data.';
end









function meth = import_illumina_infinium(filename, dataset)

global organism;

fid = fopen(filename);
data = textscan(fid, '%s %f %s %*s %*s', 'Headerlines', 1, ...
	'Delimiter', '\t', 'ReturnOnError', 0, 'TreatAsEmpty', 'NA');
fclose(fid);

sample = data{1};
beta = data{2};
gene = data{3};

samples = cell(length(unique(sample)), 1);

methylation = struct;
methylation.Degree = nan(length(organism.Genes.Name), length(samples));

gene_indices = gene_idx(gene);

current_sample = '';
sample_idx = 0;

unmapped_genes = 0;

run_ends = [ find(~strcmp(gene(1:end-1), gene(2:end))); length(gene) ];
run_lengths = diff([0; run_ends]);

pos = 1;
for r = 1:length(run_lengths)
	if ~strcmp(sample{pos}, current_sample)
		current_sample = sample{pos};
		sample_idx = sample_idx + 1;
		samples{sample_idx} = current_sample;
	end
	
	idx = gene_indices(pos);
	if isnan(idx)
		if sample_idx == 1
			unmapped_genes = unmapped_genes + 1;
		end
	else
		methylation.Degree(idx, sample_idx) = ...
			nanmean(beta(pos:pos+run_lengths(r)-1));
	end
	
	pos = pos + run_lengths(r);
end

fprintf(1, '%d genes were not mapped.\n', unmapped_genes);

methylation.Meta = struct;
methylation.Meta.Type = 'DNA methylation';
methylation.Meta.Sample = samples;
methylation.Meta.Platform = 'Illumina Infinium Human DNA Methylation 27';

[methylation.Meta.Patient, methylation.Meta.Misc] = ...
	tcga_patient_metadata(methylation.Meta.Sample);

create_dataset('TCGA GBM Illumina Infinium methylation', methylation);








function [] = import_illumina_oma(filename, dataset, array_ver)

global organism;

fid = fopen(sprintf(['/worktmp/annalam/' ...
	'jhu-usc.edu__IlluminaDNAMethylation_OMA00%d_CPI__beta-value.txt'], ...
	array_ver));

header = fgetl(fid);
if header == -1, return, end

column_names = regexp(header, '(.+?)(\t|$)', 'tokens');
num_samples = length(column_names) - 1

samples = cell(num_samples, 1);
for k = 1:num_samples
	column_name = column_names{k + 1};
	samples{k} = column_name{1};
end

data = textscan(fid, ['%s' repmat('%f', 1, num_samples)], 'Headerlines', 2, ...
	'Delimiter', '\t', 'ReturnOnError', 0, 'TreatAsEmpty', 'N/A');

gene = data{1};
if array_ver ~= 3
	for k = 1:length(gene)
		comp_name = gene{k};
		pos = find(comp_name == '_');
		if isempty(pos), continue, end
		gene{k} = comp_name(1:pos-1);
	end
elseif array_ver == 3
	fid = fopen(['/worktmp/annalam/' ...
		'jhu-usc.edu_TCGA_IlluminaDNAMethylation_OMA003_CPI.adf.txt']);
	tmp = textscan(fid, '%s %*s %s %*s %*s %*s', 'Headerlines', 1, ...
		'Delimiter', '\t');
	fclose(fid);

	probe_map = containers.Map(tmp{1}, tmp{2});
	gene = probe_map.values(gene);
end

methylation = struct;
methylation.Degree = nan(length(organism.Genes.Name), num_samples);

gene_indices = gene_idx(gene);
unmapped_genes = sum(isnan(gene_indices));

for s = 1:num_samples
	beta = data{s + 1};
	for k = 1:length(gene_indices)
		idx = gene_indices(k);
		if ~isnan(idx)
			methylation.Degree(idx, s) = beta(k);
		end
	end
end

fprintf(1, '%d genes were not mapped.\n', unmapped_genes);

methylation.Meta = struct;
methylation.Meta.Type = 'DNA methylation';
methylation.Meta.Sample = samples;
methylation.Meta.Platform = repmat( ...
	{sprintf('Illumina DNA Methylation OMA00%d Cancer Panel I', array_ver)}, ...
	num_samples, 1);
	
[methylation.Meta.Patient, methylation.Meta.Misc] = ...
	tcga_patient_metadata(methylation.Meta.Sample);

create_dataset(sprintf('TCGA GBM Illumina OMA %d methylation', array_ver), ...
	methylation);

	
	
	








function [] = import_tcga_mutations_maf(filename, dataset)

global organism;

fid = fopen(filename);

header = fgetl(fid);
if header == -1, error 'Missing MAF header.'; end

column_names = textscan(header, '%s', -1, 'Delimiter', '\t');
column_names = column_names{1};

for k = 1:length(column_names)
	col_name = column_names{k};
	if strcmpi(col_name, 'Entrez_Gene_Id'), entrez_col = k; end
	if strcmpi(col_name, 'Start_Position'), start_pos_col = k; end
	%if strcmpi(col_name, 'End_Position'), end_pos_col = k; end
	if strcmpi(col_name, 'Chrom'), chr_col = k; end
	if strcmpi(col_name, 'Strand'), strand_col = k; end
	if strcmpi(col_name, 'Center'), center_col = k; end
	if strcmpi(col_name, 'Variant_Classification'), class_col = k; end
	if strcmpi(col_name, 'Variant_Type'), type_col = k; end
	if strcmpi(col_name, 'Tumor_Sample_Barcode'), sample_col = k; end
	if strcmpi(col_name, 'Match_Norm_Sample_Barcode'), normal_col = k; end
	if strcmpi(col_name, 'Validation_Status'), validation_col = k; end
	if strcmpi(col_name, 'Mutation_Status'), mutation_status_col = k; end
	if strcmpi(col_name, 'Reference_Allele'), ref_allele_col = k; end
	if strcmpi(col_name, 'Tumor_Seq_Allele1'), mut_allele1_col = k; end
	if strcmpi(col_name, 'Tumor_Seq_Allele2'), mut_allele2_col = k; end
end

% Construct a format string for parsing.
parse_format = repmat('%s', 1, length(column_names));
data = textscan(fid, parse_format, 'Delimiter', '\t', 'BufSize', 16384);

entrez = str2double(data{entrez_col});
start_pos = str2double(data{start_pos_col});
%end_pos = str2double(data{end_pos_col});
chr = chromosome_sym2num(data{chr_col});
strand = data{strand_col};
source = data{center_col};
variant_class = data{class_col};
variant_type = data{type_col};
sample = data{sample_col};
normal = data{normal_col};
validation = data{validation_col};
mutation_status = data{mutation_status_col};
ref_allele = data{ref_allele_col};
mut_allele = cat(2, data{mut_allele1_col}, data{mut_allele2_col});

% Use liftOver to convert from hg18 coordinates to hg19 coordinates.
hg18_file = ptemp;
hg19_file = ptemp;
unmapped_file = ptemp;
fid = fopen(hg18_file, 'W');
for k = 1:length(chr)
	fprintf(fid, '%s\t%d\t%d\n', ['chr' organism.Chromosomes.Name{chr(k)}], ...
		start_pos(k) - 1, start_pos(k));
end
fclose(fid);

status = unix(sprintf(['%s/tools/liftover/liftOver %s ' ...
	'%s/tools/liftover/hg18ToHg19.over.chain %s %s'], ...
	ppath, hg18_file, ppath, hg19_file, unmapped_file));
if status ~= 0, error 'Liftover failed.'; end

fid = fopen(hg19_file);
data = textscan(fid, '%s %d %d');
fclose(fid);

if length(data{1}) ~= length(chr)
	error('Liftover only mapped %d out of %d coordinates.', length(data{1}), ...
		length(chr));
end

chr = chromosome_sym2num(data{1});
start_pos = data{2} + 1;

safe_delete(hg18_file);
safe_delete(hg19_file);
safe_delete(unmapped_file);

% Separate mutations by sample.
entrez_to_gene = containers.Map(organism.Genes.EntrezID, organism.Genes.Name);
sample_to_idx = containers.Map;
sample_count = 0;

for k = 1:length(sample)
	sample{k} = sample{k}(1:15);
	normal{k} = normal{k}(1:15);
end

mutations = struct;
mutations.Mutations = {};

mutations.Meta = struct;
mutations.Meta.Type = 'Mutations';
mutations.Meta.Sample = struct;
mutations.Meta.Sample.ID = {};

mutations.Meta.Ref = struct;
mutations.Meta.Ref.Sample = struct;
mutations.Meta.Ref.Sample.ID = {};

for k = 1:length(sample)
	if ~sample_to_idx.isKey(sample{k})
		sample_count = sample_count + 1;
		sample_to_idx(sample{k}) = sample_count;
		
		idx = find(strcmp(sample{k}, sample));
		
		mut = struct;
		mut.RefAllele = ref_allele(idx);
		mut.Reference = strcat('chr', organism.Chromosomes.Name(chr(idx)));
		
		valid = entrez_to_gene.isKey(num2cell(entrez(idx)));
		if any(valid)
			mut.Gene(valid, 1) = ...
				entrez_to_gene.values(num2cell(entrez(idx(valid))));
		end
		mut.Gene(~valid, 1) = repmat({'-'}, sum(~valid), 1);
		
		mut.Offset = start_pos(idx);
		mut.Validation = validation(idx);
		mut.Source = source(idx);
		
		use_mut_col = strcmp(ref_allele(idx), mut_allele(idx, 1));
		mut.MutAllele = mut_allele(idx + use_mut_col * size(mut_allele, 1));
		
		mutations.Mutations{sample_count} = mut;
		mutations.Meta.Sample.ID{sample_count, 1} = sample{k};
		mutations.Meta.Ref.Sample.ID{sample_count, 1} = normal{k};
	end
end

create_dataset(dataset, mutations);
augment_meta_tcga(dataset);




