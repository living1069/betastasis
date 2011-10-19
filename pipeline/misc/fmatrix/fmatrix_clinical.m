function fmatrix = fmatrix_clinical(dataset)

S = length(dataset.Meta.Sample.ID);

fmatrix = struct;
fmatrix.Features = {};
fmatrix.Samples = dataset.Meta.Sample.ID;
fmatrix.Data = zeros(0, S);

meta = dataset;
if isfield(dataset, 'Meta'), meta = dataset.Meta; end

if isfield(meta, 'Sample')
	if isfield(meta.Sample, 'Type')
		[fmatrix.Features{end+1, 1}, fmatrix.Data(end+1, :)] = ...
			categorical('Sample_type', meta.Sample.Type);
	end

	if isfield(meta.Sample, 'Timepoint')
		fmatrix.Features{end+1, 1} = 'N:CLIN:Timepoint';
		fmatrix.Data(end+1, :) = meta.Sample.Timepoint;
	end
end






function [fname, values] = categorical(name, data)

[uniq_val, ~, values] = unique(data);
fname = ['C:CLIN:' name ':'];
for k = 1:length(uniq_val)-1
	fname = [fname uniq_val{k} ','];
end
fname = [fname uniq_val{end}];

