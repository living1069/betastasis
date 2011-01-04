
% UARRAY_EXPRESSION_RMA  Calculate expression values for transcriptomic features
%
%    EXPR = UARRAY_EXPRESSION_RMA(SAMPLES, PROBESETS) calculates summarized
%    expression values for the transcriptomic features probed by PROBESETS.
%    The expression values are calculated using RMA based on the microarray 
%    probe intensities in SAMPLES.
%
%    All expression values are returned in the natural scale and are not log-
%    transformed. The rows of the expression matrix in EXPR correspond
%    one-to-one with the probesets in PROBESETS. The probesets also determine
%    the type of the expression values returned.
%
%    Probesets accepted by this  function include:
%    - Gene expression probesets
%    - Transcript expression probesets
%    - Exon expression probesets
%    - MicroRNA probesets

% Author: Matti Annala <matti.annala@tut.fi>

function expr = uarray_expression_rma(samples, probesets)

global organism;

S = size(samples.Mean, 2);

if ~iscellstr(samples.Meta.Platform)
	samples.Meta.Platform = repmat({samples.Meta.Platform}, S, 1);
end

platform = samples.Meta.Platform{1};

	
	
	
	
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
	
	if isfield(samples.Meta, 'Platform')
		expr.Meta.Platform = samples.Meta.Platform;
	end
end

