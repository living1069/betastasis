function fmatrix = cat_fmatrix(varargin)

fmatrix.data = [];
fmatrix.features = {};

matrices = varargin;

samples = {};
for k = 1:length(matrices)
	samples = cat(1, samples, matrices{k}.samples);
end

samples = unique_preserve_order(samples);

S = length(samples);
D = length(matrices);

fmatrix.samples = samples;
fmatrix.data = nan(0, S);

% samples_found{d} represents the samples of dataset d that are found in the
% merged sample list. sample_indices{d} specifies the indices in dataset d
% where those samples can be found.
samples_found = cell(D, 1);
sample_indices = cell(D, 1);

F = 0;

for d = 1:D
	m = matrices{d};
	sample_map = containers.Map(m.samples, 1:length(m.samples));
	
	samples_found = sample_map.isKey(samples);
	sample_indices = cell2mat(sample_map.values(samples(samples_found)));
	
	fmatrix.features = cat(1, fmatrix.features, m.features);
	fmatrix.data(F+1:F+length(m.features), :) = nan(length(m.features), S);
	fmatrix.data(F+1:F+length(m.features), samples_found) = ...
		m.data(:, sample_indices);
		
	F = F + length(m.features);
end

% FIXME: Add support for merging replicated feature fields.

