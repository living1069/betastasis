function [] = export_betastasis(datasets, path)

global organism;

if isstruct(datasets), datasets = { datasets }; end
	
gene_expr_dataset = [];
mirna_expr_dataset = [];
cna_seg_dataset = [];

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






	


if ~isempty(gene_expr_dataset)
	ds = datasets{gene_expr_dataset};
	
	genes = find(~any(isnan(ds.Mean), 2));
	
	log_expr = zeros(length(genes), S);
	log_expr(:, samples_found{gene_expr_dataset}) = log2(ds.Mean(genes, ...
		sample_indices{gene_expr_dataset}));
		
	for g = 1:length(genes)
		idx = genes(g);
		name = organism.Genes.Name{idx};
		dir = [path '/gene_expr/' lower(name(1))];
		[~, ~] = mkdir(dir);
		
		export_json([dir '/' name '.json'], 'expr', log_expr(g, :));
	end

	export_json([path '/gene_expr/genes.json'], ...
		'genes', organism.Genes.Name(genes));
end





if ~isempty(mirna_expr_dataset)
	ds = datasets{mirna_expr_dataset};
	
	mirnas = find(~any(isnan(ds.Mean), 2));
	
	log_expr = zeros(length(mirnas), S);
	log_expr(:, samples_found{mirna_expr_dataset}) = log2(ds.Mean(mirnas, ...
		sample_indices{mirna_expr_dataset}));
	
	dir = [path '/mirna_expr'];
	[~, ~] = mkdir(dir);

	for g = 1:length(mirnas)
		idx = mirnas(g);
		name = organism.miRNA.Name{idx};
		
		export_json([dir '/' name '.json'], 'expr', log_expr(g, :));
	end
	
	export_json([path '/mirna_expr/mirnas.json'], ...
		'mirnas', organism.miRNA.Name(mirnas));
end







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

[~, ~] = mkdir([path '/clinical']);

export_json([path '/clinical/sample_id.json'], 'values', sample_id);
export_json([path '/clinical/sample_type.json'], 'values', sample_type);

export_json([path '/clinical/patient_id.json'], 'values', patient_id);
export_json([path '/clinical/survival.json'], ...
	'survival', int32(survival), ...
	'censored', logical(censored));
export_json([path '/clinical/total_gleason.json'], ...
	'total_gleason', int32(total_gleason));

	

	
	
	
	
	
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

