function fmatrix = fmatrix_clinical(dataset)

meta = dataset;
if isfield(dataset, 'meta'), meta = dataset.meta; end

S = length(meta.sample_id);

fmatrix.features = {};
fmatrix.samples = meta.sample_id;
fmatrix.data = zeros(0, S);

if isfield(meta, 'sample_type')
	[fmatrix.features{end+1, 1}, fmatrix.data(end+1, :)] = ...
		categorical('sample_type', meta.sample_type);
end

if isfield(meta, 'survival_time')
	fmatrix.features{end+1, 1} = 'N:CLIN:survival_time';
	fmatrix.data(end+1, :) = meta.survival_time;
end

if isfield(meta, 'survival_time_censored')
	fmatrix.features{end+1, 1} = 'N:CLIN:survival_censored';
	fmatrix.data(end+1, :) = meta.survival_time_censored;
end

if isfield(meta, 'progression_time')
	fmatrix.features{end+1, 1} = 'N:CLIN:progression_time';
	fmatrix.data(end+1, :) = meta.progression_time;
end

if isfield(meta, 'progression_censored')
	fmatrix.features{end+1, 1} = 'N:CLIN:progression_censored';
	fmatrix.data(end+1, :) = meta.progression_censored;
end

if isfield(meta, 'timepoint')
	[fmatrix.features{end+1, 1}, fmatrix.data(end+1, :)] = ...
		categorical('timepoint', meta.timepoint);
end






function [fname, values] = categorical(name, data)

[uniq_val, ~, values] = unique(data)
fname = ['C:CLIN:' name ':'];
for k = 1:length(uniq_val)-1
	fname = [fname uniq_val{k} ';'];
end
fname = [fname uniq_val{end}];

