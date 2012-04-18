function mutations = import_mutations_maf(files)

global organism;

if nargin < 1
	files = find_files('\.maf');
end

mutations = struct;
mutations.meta.type = 'Mutations';
mutations.meta.sample_id = {};
mutations.chromosome = {};
mutations.position = {};
mutations.ref_allele = {};
mutations.mut_allele = {};
mutations.genotype = [];

for f = 1:length(files)
	[data, headers] = readtable(files{f});

	start_pos = str2double(data{rx(headers, 'start.*pos')});
	chr = chromosome_sym2num(data{rx(headers, 'chrom')});
	strand = data{rx(headers, 'strand')};
	sample = data{rx(headers, 'tumor.*sample.*barcode')};
	normal = data{rx(headers, 'norm.*sample.*barcode')};
	ref_allele = data{rx(headers, 'ref.*allele')};
	mut_allele1 = data{rx(headers, 'tumor.*seq.*allele1')};
	mut_allele2 = data{rx(headers, 'tumor.*seq.*allele2')};

	% Use liftOver to convert from hg18 coordinates to hg19 coordinates.
	if do_liftover
	end

	% Separate mutations by sample.
	entrez_to_gene = containers.Map(organism.Genes.EntrezID, organism.Genes.Name);
	sample_to_idx = containers.Map;
	sample_count = 0;

	for k = 1:length(sample)
		sample{k} = sample{k}(1:15);
		normal{k} = normal{k}(1:15);
	end
	
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
			
			mut.Validation = validation(idx);
			mut.Source = source(idx);
			
			use_mut_col = strcmp(ref_allele(idx), mut_allele(idx, 1));
			mut.MutAllele = mut_allele(idx + use_mut_col * size(mut_allele, 1));
			
			mutations.Mutations{sample_count} = mut;
			mutations.Meta.Sample.ID{sample_count, 1} = sample{k};
			mutations.Meta.Ref.Sample.ID{sample_count, 1} = normal{k};
		end
	end
end



