function mutations = read_mutations_tab(filename, varargin)

global organism;
genes = organism.Genes;

method = {};

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Method')
		method = varargin{k+1};
		if iscell(method) && size(method, 2) == 2
			method = containers.Map(method(:, 1), method(:, 2));
		elseif isa(method, 'containers.Map')
			method = method;
		end
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

fid = fopen(filename);

header = fgetl(fid);
if header == -1, error 'Missing header line.'; end

column_names = textscan(header, '%s', -1, 'Delimiter', '\t');
column_names = column_names{1};

% Construct a format string for parsing.
parse_format = repmat('%s', 1, length(column_names));
data = textscan(fid, parse_format, 'Delimiter', '\t', 'BufSize', 16384);

mut = struct;

for k = 1:length(column_names)
	col_name = column_names{k};
	if regexpi(col_name, '^Entrez_Gene_Id$')
		entrez_to_gene = containers.Map(genes.EntrezID, genes.Name);

		entrez = data{k};
		valid = entrez_to_gene.isKey(entrez);
		mut.Gene = repmat({'-'}, length(valid), 1);
		mut.Gene(valid) = entrez_to_gene.values(entrez);
	end
	if regexpi(col_name, '^Gene$')
		idx = gene_idx(data{k});
		valid = ~isnan(idx);
		for k = find(~valid)'
			fprintf(1, 'Unknown gene %s.\n', genes.Name{k});
		end
		mut.Gene = repmat({'-'}, length(valid), 1);
		mut.Gene(valid) = genes.Name(idx(valid));
	end
	if regexpi(col_name, '^(Position|Start.position)$')
		mut.Position = str2double(data{k});
	end
	if regexpi(col_name, '^(Chr|Chrom|Chromosome)$')
		mut.Chromosome = chromosome_sym2num(data{k});
	end
	if regexpi(col_name, '^Strand$')
		mut.Strand = data{k};
	end
	if regexpi(col_name, '^(Tumor_Sample_Barcode|Sample)$')
		sample_ids = data{k};
	end
	if regexpi(col_name, '^Mutation.type$')
		mut.Type = data{k};
	end
	if regexpi(col_name, '^Method$')
		mut.Method = data{k};
	end
	if strcmpi(col_name, 'Center')
		mut.Source = data{k};
	end
	if strcmpi(col_name, 'Variant_Classification')
		mut.VariantClass = data{k};
	end
	if strcmpi(col_name, 'Variant_Type')
		mut.VariantType = data{k};
	end
	if strcmpi(col_name, 'Match_Norm_Sample_Barcode')
		mut.NormalSampleID = data{k};
	end
	if strcmpi(col_name, 'Validation_Status')
		mut.Validation = data{k};
	end
	if strcmpi(col_name, 'Mutation_Status')
		mut.Status = data{k};
	end
	if regexpi(col_name, '^(Reference_Allele|Reference)$')
		mut.RefAllele = data{k};
	end
	if regexpi(col_name, '^(Mutant|Mutant.?allele|Tumor_Seq_Allele1?)$')
		mut.MutAllele = data{k};
	end
end

mut = liftover(mut, 'hg18 -> hg19');

% Separate mutations by sample.
mutations = struct;
mutations.Mutations = {};

mutations.Meta = struct;
mutations.Meta.Type = 'Genomic variation';
mutations.Meta.Sample = struct;
mutations.Meta.Sample.ID = {};

mutations.Meta.Ref = struct;
mutations.Meta.Ref.Sample = struct;
mutations.Meta.Ref.Sample.ID = {};

% If we have method information available, we need to include more genes
if isempty(method)
	uniq_samples = unique(sample_ids);
else
	uniq_samples = unique(cat(1, sample_ids, method.keys'));
	mutations.Method = repmat({'None'}, 1, length(uniq_samples));
	valid = method.isKey(uniq_samples);
	mutations.Method(valid) = method.values(uniq_samples(valid));
end

for s = 1:length(uniq_samples)
	idx = find(strcmp(uniq_samples{s}, sample_ids))
	
	mutations.Mutations{1, s} = filter_struct(mut, idx);
	mutations.Meta.Sample.ID{s, 1} = uniq_samples{s};
	%mutations.Meta.Ref.Sample.ID{sample_count, 1} = normal{k};
end


