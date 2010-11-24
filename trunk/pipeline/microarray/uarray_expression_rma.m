
% This function calculates the expression of genes, transcripts or exons
% based on a set of microarray experiment replicates. The microarray
% experiments must be executed on identical platforms.
% 
% A probeset needs to be provided in order for summarization of the probe
% intensities to be possible.
%
% Inputs:
%     samples - An NxM matrix where each column represents a microarray
%         experiment and each row represents the intensities of a particular
%         probe in the set of replicate experiments.
%     probesets - A data structure that describes how probes are gathered into
%         into probesets. This data structure plays a role in summarization.
%
% Outputs:
%     expr - Microarray expression results
%
% Author: Matti Annala <matti.annala@tut.fi>

function expr = uarray_expression_rma(samples, probesets)

global organism;

platform = samples.Meta.Platform;
if iscellstr(platform), platform = platform{1}; end

	
	
	
	
if any(any(samples.Mean <= 0))
	fprintf(1, ['Data contains zero or negative intensities. ' ...
		'Replacing with the lowest positive intensity...\n']);
	samples.Mean(samples.Mean <= 0) = Inf;
	samples.Mean(samples.Mean == Inf) = min(min(samples.Mean));
end




	
	
if regexp(platform, 'Agilent')
	fprintf(1, 'Skipping background adjustment on Agilent microarray...\n');
else
	fprintf(1, 'Performing RMA background adjustment on samples...\n');
	samples.Mean = rmabackadj(samples.Mean);
end





fprintf(1, 'Normalizing probe intensities using quantile normalization...\n');

progress = Progress;

%sketch_cols = 1:floor((size(samples.Mean, 2)-1)/S):size(samples.Mean, 2);
sketch_mean = zeros(size(samples.Mean, 1), 1, 'single');

S = size(samples.Mean, 2);

for s = 1:S
	[sorted, samples.Mean(:, s)] = sort(samples.Mean(:, s), 'descend');
	sketch_mean = sketch_mean + sorted / S;
	progress.update(0.7 * s/S);
end
for s = 1:S
	samples.Mean(samples.Mean(:, s), s) = sketch_mean;
	progress.update(0.7 + 0.3 * s/S);
end





fprintf(1, 'Summarizing expression data using RMA summarization...\n');

expr = struct( ...
	'Mean', zeros(length(probesets.ProbeCount), size(samples.Mean, 2)));

% Reorder the probes so that the probes of each probeset are placed contiguously
% in the vector. The probe index vector therefore has the following form:
% 
%   [ 0 1 2 0 1 0 1 2 3 4 5 6 ]'
% 
% In the example index vector above, we have three probesets with 3, 2 and
% 7 probes, respectively.

probe_reordering = zeros(sum(probesets.ProbeCount), 1);
probe_indices = zeros(sum(probesets.ProbeCount), 1);
k = 1;
for ps = 1:length(probesets.ProbeCount)
	pc = probesets.ProbeCount(ps);
	probe_reordering(k:k+pc-1) = probesets.Probes(ps, 1:pc);
	probe_indices(k:k+pc-1) = 0:pc-1;
	k = k + pc;
end

rma_ordered_sample_data = samples.Mean(probe_reordering, :);
summary = rmasummary(probe_indices, rma_ordered_sample_data, ...
	'Output', 'natural');

% Set the expression of empty probesets to NaN.
expr.Mean(:, :) = NaN;
expr.Mean(find(probesets.ProbeCount ~= 0), :) = summary;
if isfield(samples, 'Meta')
	expr.Meta = samples.Meta;
	expr.Meta.Type = probesets.Type;
	expr.Meta.SummarizationMethod = repmat({'RMA'}, size(samples.Mean, 2), 1);
	expr.Meta.Organism = probesets.Organism;
	
	if isfield(probesets, 'GenomeVersion')
		expr.Meta.GenomeVersion = probesets.GenomeVersion;
	elseif isfield(probesets, 'miRNAVersion')
		expr.Meta.miRNAVersion = probesets.miRNAVersion;
	end
	
	expr.Meta.Platform = repmat({platform}, ...
		size(samples.Mean, 2), 1);
end

