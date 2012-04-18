function [fmatrix, gcna] = feature_matrix(varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;
exons = organism.Exons;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

fmatrix = struct;
fmatrix.Data = [];
fmatrix.Samples = {};
fmatrix.Features = {};
	
gene_expr_dataset = [];
exon_expr_dataset = [];
mirna_expr_dataset = [];
cna_dataset = [];
methylation_dataset = [];
mutation_dataset = [];

datasets = varargin;
for k = 1:length(datasets)
	if ~isfield(datasets{k}, 'Meta')
		error('Dataset #%d has no metadata.', k);
	end
	
	type = datasets{k}.Meta.Type;
	if strcmpi(type, 'Gene expression')
		if ~isempty(gene_expr_dataset)
			error 'More than one gene expression dataset found.';
		end
		gene_expr_dataset = k;
	elseif strcmpi(type, 'Exon expression')
		if ~isempty(exon_expr_dataset)
			error 'More than one exon expression dataset found.';
		end
		exon_expr_dataset = k;
	elseif strcmpi(type, 'MicroRNA expression') || ...
		strcmpi(type, 'miRNA expression')
		
		if ~isempty(mirna_expr_dataset)
			error 'More than one microRNA expression dataset found.';
		end
		mirna_expr_dataset = k;
	elseif strcmpi(type, 'Copy number segments')
		if ~isempty(cna_dataset), error 'More than one CNA dataset found.'; end
		cna_dataset = k;
	elseif strcmpi(type, 'DNA methylation') || strcmpi(type, 'Gene methylation')
		if ~isempty(methylation_dataset)
			error 'More than one methylation dataset found.';
		end
		methylation_dataset = k;
	elseif strcmpi(type, 'Genomic variation')
		if ~isempty(mutation_dataset)
			error 'More than one mutation dataset found.';
		end
		mutation_dataset = k;
	else
		error('%s datasets are not supported.', datasets{k}.Meta.Type);
	end
end

samples = [];
for k = 1:length(datasets)
	ds = datasets{k};
	ds_samples = ds.Meta.Sample.ID;
	if length(ds_samples) ~= length(unique(ds_samples))
		fprintf(1, ['Dataset #%d contains technical replicates. ' ...
			'Merging the replicates now...\n'], k);
		mds = merge_replicates(ds);
		ds_samples = ds.Meta.Sample.ID;
	end
	samples = [samples; ds_samples];
end

samples = unique(samples);
samples = samples(~strcmp('-', samples));
S = length(samples);
D = length(datasets);

% samples_found{d} represents the samples of dataset d that are found in the
% merged sample list. sample_indices{d} specifies the indices in dataset d
% where those samples can be found.
samples_found = cell(D, 1);
sample_indices = cell(D, 1);

for d = 1:D
	ds = datasets{d};
	valid = find(~strcmp('-', ds.Meta.Sample.ID));
	sample_map = containers.Map(ds.Meta.Sample.ID(valid), num2cell(valid));
	
	samples_found{d} = sample_map.isKey(samples); 
	sample_indices{d} = cell2mat(sample_map.values(samples(samples_found{d})));
end

fmatrix.Samples = samples;
F = 0;





	


if ~isempty(gene_expr_dataset)
	ds = datasets{gene_expr_dataset};
	
	valid = find(~any(isnan(ds.Mean), 2));
	
	fmatrix.Features = cat(1, fmatrix.Features, cell(length(valid), 1));
	for g = 1:length(valid)
		if isnan(genes.Chromosome(valid(g)))
			pos_str = '?:?:?:?';
		else
			pos_str = sprintf('%s:%d:%d:%s', ...
				chromosomes.Name{genes.Chromosome(valid(g))}, ...
				genes.Position(valid(g), 1), genes.Position(valid(g), 2), ...
				genes.Strand(valid(g)));
		end
		
		fmatrix.Features{F+g} = ...
			sprintf('N:EXPR:%s:%s', genes.Name{valid(g)}, pos_str);
	end
	
	log_expr = log2(ds.Mean(valid, sample_indices{gene_expr_dataset}));
	if any(any(log_expr == -Inf))
		fprintf(1, ['WARNING: Negative infinities found in gene expression '...
			'data. Replacing with lowest finite value...\n']);
		log_expr(log_expr == -Inf) = Inf;
		log_expr(log_expr == Inf) = nanmin(nanmin(log_expr));
	end
	
	fmatrix.Data = [fmatrix.Data; nan(length(valid), S)];
	fmatrix.Data(F+1:F+length(valid), samples_found{gene_expr_dataset}) = ...
		log_expr;
	
	fprintf(1, '%d gene expression features.\n', length(valid));
	F = F + length(valid);
end








if ~isempty(gene_expr_dataset) && ~isempty(exon_expr_dataset)
	gene_ds = datasets{gene_expr_dataset};
	exon_ds = datasets{exon_expr_dataset};
	
	% We must pair the gene and exon expression samples in order to calculate
	% gene-normalized exon expression values.
	gene_log_expr = nan(size(gene_ds.Mean, 1), S);
	exon_log_expr = nan(size(exon_ds.Mean, 1), S);
	
	gene_log_expr(:, samples_found{gene_expr_dataset}) = ...
		log2(gene_ds.Mean(:, sample_indices{gene_expr_dataset}));
	exon_log_expr(:, samples_found{exon_expr_dataset}) = ...
		log2(exon_ds.Mean(:, sample_indices{exon_expr_dataset}));
		
	splice_index = nan(size(exon_log_expr));
	
	fprintf(1, 'WARNING: Cheap trick!\n');
	
	for e = 1:size(exon_log_expr, 1)
		g = exons.Gene(e);
		splice_index(e, :) = exon_log_expr(e, :) - gene_log_expr(g, :);
		splice_index(e, exon_log_expr(e, :) < 3) = NaN;  % FIXME
	end
	
	%valid = find(any(~isnan(splice_index), 2));
	
	E = size(exon_log_expr, 1);
	%fmatrix.Features = cat(1, fmatrix.Features, cell(length(valid), 1));
	fmatrix.Features = cat(1, fmatrix.Features, cell(E, 1));
	%for f = 1:length(valid)
	for ex = 1:E
		%ex = valid(f);
		g = exons.Gene(ex);
		if isnan(genes.Chromosome(g))
			pos_str = '?:?:?:?';
		else
			pos_str = sprintf('%s:%d:%d:%s', ...
				chromosomes.Name{genes.Chromosome(g)}, ...
				genes.Position(g, 1), genes.Position(g, 2), genes.Strand(g));
		end
		
		fmatrix.Features{F+ex} = sprintf('N:ESPL:%s[%s]:%s', genes.Name{g}, ...
			exons.ID{ex}, pos_str);
	end
	
	fmatrix.Data = [fmatrix.Data; nan(E, S)];
	fmatrix.Data(F+1:F+E, :) = splice_index;
	%fmatrix.Data(F+1:F+length(valid), :) = splice_index(valid, :);
	
	fprintf(1, '%d exon splicing features.\n', E);
	F = F + E;
end








if ~isempty(mirna_expr_dataset)
	ds = datasets{mirna_expr_dataset};
	
	if isfield(ds, 'Mean')
		expr = ds.Mean;
	elseif isfield(ds, 'miRNA')
		expr = ds.miRNA;
	end
	
	valid = find(~any(isnan(expr), 2));
	
	fmatrix.Features = cat(1, fmatrix.Features, cell(length(valid), 1));
	for m = 1:length(valid)
		[pre, ~] = find(pre_mirnas.Matures == valid(m));
		if isempty(pre) || isnan(pre_mirnas.Chromosome(pre(1)))
			pos_str = '?:?:?:?';
		else
			pre = pre(1);
			pos_str = sprintf('%s:%d:%d:%s', ...
				chromosomes.Name{pre_mirnas.Chromosome(pre)}, ...
				pre_mirnas.Position(pre, 1), pre_mirnas.Position(pre, 2), ...
				pre_mirnas.Strand(pre));
		end

		fmatrix.Features{F+m} = sprintf('N:EXPR:%s:%s', ...
			mirnas.Name{valid(m)}, pos_str);
	end
	
	fmatrix.Data = [fmatrix.Data; nan(length(valid), S)];
	fmatrix.Data(F+1:F+length(valid), samples_found{mirna_expr_dataset}) = ...
		log2(expr(valid, sample_indices{mirna_expr_dataset}));
	
	fprintf(1, '%d microRNA expression features.\n', length(valid));
	F = F + length(valid);
end






if ~isempty(cna_dataset)
	ds = datasets{cna_dataset};
	
	gcna = gene_cna(ds, 'MinRegion', 100e3);
	
	valid = find(~any(isnan(gcna), 2));
	
	fmatrix.Features = cat(1, fmatrix.Features, cell(length(valid), 1));
	for g = 1:length(valid)
		if isnan(genes.Chromosome(valid(g)))
			pos_str = '?:?:?:?';
		else
			pos_str = sprintf('%s:%d:%d:%s', ...
				chromosomes.Name{genes.Chromosome(valid(g))}, ...
				genes.Position(valid(g), 1), genes.Position(valid(g), 2), ...
				genes.Strand(valid(g)));
		end

		fmatrix.Features{F+g} = ...
			sprintf('N:CNA:%s:%s', genes.Name{valid(g)}, pos_str);
	end
	
	fmatrix.Data = [fmatrix.Data; nan(length(valid), S)];
	fmatrix.Data(F+1:F+length(valid), samples_found{cna_dataset}) = ...
		gcna(valid, sample_indices{cna_dataset});
	
	fprintf(1, '%d copy number features.\n', length(valid));
	F = F + length(valid);
end







if ~isempty(methylation_dataset)
	ds = datasets{methylation_dataset};
	
	valid = find(~any(isnan(ds.Degree), 2));
	
	fmatrix.Features = cat(1, fmatrix.Features, cell(length(valid), 1));
	for g = 1:length(valid)
		if isnan(genes.Chromosome(valid(g)))
			pos_str = '?:?:?:?';
		else
			pos_str = sprintf('%s:%d:%d:%s', ...
				chromosomes.Name{genes.Chromosome(valid(g))}, ...
				genes.Position(valid(g), 1), genes.Position(valid(g), 2), ...
				genes.Strand(valid(g)));
		end

		fmatrix.Features{F+g} = ...
			sprintf('N:METH:%s:%s', genes.Name{valid(g)}, pos_str);
	end
	
	fmatrix.Data = [fmatrix.Data; nan(length(valid), S)];
	fmatrix.Data(F+1:F+length(valid), samples_found{methylation_dataset}) = ...
		ds.Degree(valid, sample_indices{methylation_dataset});
	
	fprintf(1, '%d DNA methylation features.\n', length(valid));
	F = F + length(valid);
end








if ~isempty(mutation_dataset)
	ds = datasets{mutation_dataset};
	
	all_mut_strs = {};
	mut_strs = {};
	for s = 1:length(ds.Mutations)
		mut_strs{s} = {};
		for m = 1:length(ds.Mutations{s}.Gene)
			mut = ds.Mutations{s};
			mut_strs{s}{m} = sprintf('N:MUT:%s:%s:%d:%d:?', mut.Gene{m}, ...
				chromosomes.Name{mut.Chromosome(m)}, mut.Position(m), ...
				mut.Position(m));
			all_mut_strs{end+1, 1} = mut_strs{s}{m};
		end
	end
	
	all_mut_strs = unique(all_mut_strs);
	
	mut_matrix = zeros(length(all_mut_strs), length(mut_strs));
	for s = 1:length(mut_strs)
		[~, loc] = ismember(mut_strs{s}, all_mut_strs);
		mut_matrix(loc(loc ~= 0), s) = 1;
	end
	
	% Some samples may have no mutation information. Make sure that we place
	% only NaN values in those columns.
	mut_matrix(:, strcmpi('None', ds.Method)) = NaN;
	
	fmatrix.Features(F+1:F+length(all_mut_strs), 1) = all_mut_strs;
	
	fmatrix.Data = [fmatrix.Data; nan(length(all_mut_strs), S)];
	fmatrix.Data(F+1:F+length(all_mut_strs), ...
		samples_found{mutation_dataset}) = ...
		mut_matrix(:, sample_indices{mutation_dataset});
		
	fprintf(1, '%d mutation features.\n', length(all_mut_strs));
	F = F + length(all_mut_strs);
end








gender = nan(S, D);
survival = nan(S, D);
censored = nan(S, D);
rft = nan(S, D);
rft_censored = nan(S, D);

for d = 1:D
	meta = datasets{d}.Meta;
	
	found = samples_found{d};
	idx = sample_indices{d};
		
	if isfield(meta, 'Patient')
		if isfield(meta.Patient, 'Gender')
			male = strcmpi('Male', meta.Patient.Gender(idx));
			female = strcmpi('Female', meta.Patient.Gender(idx));
			tmp = nan(length(idx), 1); tmp(male) = 1; tmp(female) = 0;
			gender(found, d) = tmp;
		end
		
		if isfield(meta.Patient, 'SurvivalTime')
			survival(found, d) = meta.Patient.SurvivalTime(idx);
			censored(found, d) = meta.Patient.Censored(idx);
		end
		
		if isfield(meta.Patient, 'RecurrenceFreeTime')
			rft(found, d) = meta.Patient.RecurrenceFreeTime(idx);
			rft_censored(found, d) = ~meta.Patient.RecurrenceEvent(idx);
		end
	end
end

fmatrix.Features = cat(1, fmatrix.Features, ...
	{ 'C:CLIN:Gender'; ...
	  'N:CLIN:Survival time'; 'N:CLIN:Survival censored'; ...
	  'N:CLIN:Recurrence-free time'; 'N:CLIN:Recurrence censored' });
fmatrix.Data(F+5, :) = 0;

fmatrix.Data(F+1, :) = consensus('Gender', gender, samples);
fmatrix.Data(F+2, :) = consensus('Survival time', survival, samples);
fmatrix.Data(F+3, :) = consensus('Survival censored', censored, samples);
fmatrix.Data(F+4, :) = consensus('Recurrence-free time', rft, samples);
fmatrix.Data(F+5, :) = consensus('Recurrence censored', rft_censored, samples);
F = F + 5;



	

	
	
	
	
	
function ret = consensus(name, data, samples)

S = length(samples);

if isnumeric(data)
	ret = nan(size(data, 1), 1);
	for s = 1:S
		vals = unique(data(s, ~isnan(data(s, :))));
		if length(vals) > 1
			error('No consensus on %s for sample %s.', name, samples{s});
		elseif length(vals) == 1
			ret(s) = vals(1);
		end
	end
	
else
	error 'Unsupported data type.';
end	

