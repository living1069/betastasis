function [] = export_rf_matrix(datasets, path)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

if isstruct(datasets), datasets = { datasets }; end
	
gene_expr_dataset = [];
mirna_expr_dataset = [];
methylation_dataset = [];

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
	elseif strcmpi(type, 'MicroRNA expression') || ...
		strcmpi(type, 'miRNA expression')
		
		if ~isempty(mirna_expr_dataset)
			error 'More than one microRNA expression dataset found.';
		end
		mirna_expr_dataset = k;
	elseif strcmpi(type, 'DNA methylation') || strcmpi(type, 'Gene methylation')
		if ~isempty(methylation_dataset)
			error 'More than one methylation dataset found.';
		end
		methylation_dataset = k;
		
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

fprintf(1, '%d total samples.\n', S);

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





features = cell(1, 1e6);
out = nan(S, 1e6);
F = 0;   % Number of features.





if ~isempty(gene_expr_dataset)
	ds = datasets{gene_expr_dataset};
	
	valid = find(~any(isnan(ds.Mean), 2));
	
	log_expr = nan(S, length(valid));
	log_expr(samples_found{gene_expr_dataset}, :) = log2(ds.Mean(valid, ...
		sample_indices{gene_expr_dataset}))';
		
	out(:, F+1:F+length(valid)) = log_expr;
	
	for g = 1:length(valid)
		idx = valid(g);
		
		chr = 'NA';
		if ~isnan(genes.Chromosome(idx))
			chr = chromosomes.Name{genes.Chromosome(idx)};
		end
		
		features{F+g} = sprintf(['CHR=%s,DATASUPPORT=[-Inf Inf],' ...		
			'DATATYPE=CONTINUOUS,FEATURENAME=%s,FEATURETYPE=GENE,' ...
			'START=%d,STOP=%d'], chr, genes.Name{idx}, ...
			genes.Position(idx, 1), genes.Position(idx, 2));
	end
	
	fprintf(1, '%d gene expression features.\n', length(valid));
	F = F + length(valid);
end




if ~isempty(mirna_expr_dataset)
	ds = datasets{mirna_expr_dataset};
	
	valid = find(~any(isnan(ds.Mean), 2));
	
	log_expr = nan(S, length(valid));
	log_expr(samples_found{mirna_expr_dataset}, :) = log2(ds.Mean(valid, ...
		sample_indices{mirna_expr_dataset}))';
		
	out(:, F+1:F+length(valid)) = log_expr;
	
	for g = 1:length(valid)
		idx = valid(g);
		
		chr = 'NA';
		pos = NaN;
		
		[premirs, ~] = find(pre_mirnas.Matures == idx);
		premir_entrez = unique(pre_mirnas.EntrezGene(premirs));
		premir_entrez = premir_entrez(~isnan(premir_entrez));
		premir_genes = find(ismember(genes.EntrezID, premir_entrez)); 
		
		premir_chrs = unique(genes.Chromosome(premir_genes));
		premir_chrs = premir_chrs(~isnan(premir_chrs));
		
		premir_pos = unique(genes.Position(premir_genes, 1));
		premir_pos = premir_pos(~isnan(premir_pos));
		
		if length(premir_chrs) == 1
			chr = chromosomes.Name{premir_chrs};
			pos = premir_pos;
		end
		
		features{F+g} = sprintf(['CHR=%s,DATASUPPORT=[-Inf Inf],' ...		
			'DATATYPE=CONTINUOUS,FEATURENAME=%s,FEATURETYPE=MIRNA,' ...
			'START=%d,STOP=%d'], chr, mirnas.Name{idx}, pos, pos);
	end
	
	fprintf(1, '%d miRNA expression features.\n', length(valid));
	F = F + length(valid);
end







if ~isempty(methylation_dataset)
	ds = datasets{methylation_dataset};
	
	valid = find(~any(isnan(ds.Degree), 2));
	
	degree = nan(S, length(valid));
	degree(samples_found{methylation_dataset}, :) = ds.Degree(valid, ...
		sample_indices{methylation_dataset})';
		
	out(:, F+1:F+length(valid)) = degree;
	
	for g = 1:length(valid)
		idx = valid(g);
		
		chr = 'NA';
		if ~isnan(genes.Chromosome(idx))
			chr = chromosomes.Name{genes.Chromosome(idx)};
		end
		
		pos = genes.Position(idx, 1);
		
		features{F+g} = sprintf(['CHR=%s,DATASUPPORT=[0 1],' ...		
			'DATATYPE=CONTINUOUS,FEATURENAME=%s,FEATURETYPE=METH,' ...
			'START=%d,STOP=%d'], chr, genes.Name{idx}, pos, pos);
	end
	
	fprintf(1, '%d DNA methylation features.\n', length(valid));
	F = F + length(valid);
end









if 0



patient_id = repmat({'-'}, S, D);
gender = repmat({'-'}, S, D);
race = repmat({'-'}, S, D);
survival = nan(S, D);
censored = nan(S, D);
total_gleason = nan(S, D);

sample_id = repmat({'-'}, S, D);
sample_type = repmat({'-'}, S, D);

for d = 1:D
	meta = datasets{d}.Meta;
	
	found = samples_found{d};
	idx = sample_indices{d};
		
	if isfield(meta, 'Patient')
		if isfield(meta.Patient, 'ID')
			patient_id(found, d) = meta.Patient.ID(idx);
		end
		
		if isfield(meta.Patient, 'Gender')
			gender(found, d) = meta.Patient.Gender(idx);
		end
		
		if isfield(meta.Patient, 'Race')
			race(found, d) = meta.Patient.Race(idx);
		end
		
		if isfield(meta.Patient, 'SurvivalTime')
			survival(found, d) = meta.Patient.SurvivalTime(idx);
			censored(found, d) = meta.Patient.Censored(idx);
		end
	
		if isfield(meta.Patient, 'TotalGleason')
			total_gleason(found, d) = meta.Patient.TotalGleason(idx);
		elseif isfield(meta.Patient, 'PathologicalTotalGleason')
			total_gleason(found, d) = ...
				meta.Patient.PathologicalTotalGleason(idx);
		elseif isfield(meta.Patient, 'BiopsyTotalGleason')
			total_gleason(found, d) = ...
				meta.Patient.BiopsyTotalGleason(idx);
		end
	end
	
	if isfield(meta, 'Sample')
		if isfield(meta.Sample, 'ID')
			sample_id(found, d) = meta.Sample.ID(idx);
		end
		
		if isfield(meta.Sample, 'Type')
			sample_type(found, d) = meta.Sample.Type(idx);
		end
	end
end

sample_id = consensus('sample id', sample_id, samples);
sample_type = consensus('sample type', sample_type, samples);

patient_id = consensus('patient id', patient_id, samples);
survival = consensus('survival time', survival, samples);
censored = consensus('survival censoring', censored, samples);
total_gleason = consensus('total Gleason score', total_gleason, samples);

end




out_tmp = ptemp;

out = out(:, 1:F);
features = features(1:F);

dlmwrite(out_tmp, out, '\t');

tmp_fid = fopen(out_tmp);
final_fid = fopen(path, 'W');

for f = 1:length(features)
	fprintf(final_fid, '\t%s', features{f});
end
fprintf(final_fid, '\n');

for s = 1:S
	line = fgetl(tmp_fid);
	if line == -1, error 'Data file ended abruptly.', end
	
	fprintf(final_fid, '%s\t%s\n', samples{s}, line);
end

fclose(final_fid);
fclose(tmp_fid);

safe_delete(out_tmp);


	

	
	
	
	
	
function ret = consensus(name, data, samples)

S = length(samples);

if isnumeric(data)
	ret = -ones(size(data, 1), 1);
	for s = 1:S
		vals = unique(data(s, ~isnan(data(s, :))));
		if length(vals) > 1
			fprintf(1, ['No consensus on %s for sample %s. ' ...
				'Using the median value %.2f...\n'], name, samples{s}, ret(s));
			ret(s) = median(data(s, ~isnan(data(s, :))));
		elseif length(vals) == 1
			ret(s) = vals(1);
		end
	end
	
elseif iscellstr(data)
	ret = repmat({'-'}, size(data, 1), 1);
	for s = 1:S
		vals = unique(data(s, ~strcmp('-', data(s, :))));
		if length(vals) > 1
			fprintf(1, ['No consensus on %s for sample %s. ' ...
				'Leaving it empty...\n'], name, samples{s});
		elseif length(vals) == 1
			ret{s} = vals{1};
		end
	end
	
else
	error 'Unsupported data type.';
end	

