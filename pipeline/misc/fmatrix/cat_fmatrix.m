function fmatrix = cat_fmatrix(varargin)

fmatrix = struct;
fmatrix.Data = [];
fmatrix.Features = {};

matrices = varargin;

samples = {};
for k = 1:length(matrices)
	samples = cat(1, samples, matrices{k}.Samples);
end

samples = unique(samples);
S = length(samples);
D = length(matrices);

fmatrix.Samples = samples;
fmatrix.Data = nan(0, S);

% samples_found{d} represents the samples of dataset d that are found in the
% merged sample list. sample_indices{d} specifies the indices in dataset d
% where those samples can be found.
samples_found = cell(D, 1);
sample_indices = cell(D, 1);

F = 0;

for d = 1:D
	m = matrices{d};
	sample_map = containers.Map(m.Samples, 1:length(m.Samples));
	
	samples_found = sample_map.isKey(samples);
	sample_indices = cell2mat(sample_map.values(samples(samples_found)));
	
	fmatrix.Features = cat(1, fmatrix.Features, m.Features);
	fmatrix.Data(F+1:F+length(m.Features), :) = nan(length(m.Features), S);
	fmatrix.Data(F+1:F+length(m.Features), samples_found) = ...
		m.Data(:, sample_indices);
		
	F = F + length(m.Features);
end

