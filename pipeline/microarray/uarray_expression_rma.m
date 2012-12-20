
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

function expr = uarray_expression_rma(samples, probesets, varargin)

global organism;

background_adjustment = false;
quantile_normalize = true;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'background.*adj')
		background_adjustment = varargin{k+1};
		continue;
	end
	
	if rx(varargin{k}, 'normalize')
		quantile_normalize = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if ~background_adjustment && isfield(samples.meta, 'platform') && ...
	any(rx(samples.meta.platform, 'affy'))
	fprintf(['WARNING: Affymetrix array detected but no background ' ...
		'adjustment requested.\n']);
end

S = size(samples.mean, 2);


	
	
if any(any(samples.mean <= 0))
	fprintf(['Data contains zero or negative intensities. ' ...
		'Replacing with the lowest positive intensity...\n']);
	samples.mean(samples.mean <= 0) = Inf;
	samples.mean(samples.mean == Inf) = min(min(samples.mean));
end

if any(any(isnan(samples.mean))) || any(any(samples.mean == Inf))
	error 'Probe intensity data contains NaN or Inf values.';
end


	
	
if background_adjustment
	fprintf('Performing RMA background adjustment on samples...\n');
	samples.mean = rmabackadj(samples.mean);
	
	% Sometimes RMA background adjustment introduces Inf or NaN values in data.
	% Check for such samples, and warn the user about them.
	bad = find(any(isnan(samples.mean), 1) | any(samples.mean == Inf, 1));
	if length(bad) > 0
		bad_str = sprintf('%d', bad(1));
		for k = 2:length(bad)
			bad_str = sprintf('%s, %d', bad_str, bad(k));
		end
		error(['Background adjustment seems to have introduced Inf or NaN ' ...
			'values in samples %s.'], bad_str);
	end
end




if quantile_normalize
	fprintf(1, ...
		'Normalizing probe intensities using quantile normalization...\n');

	progress = Progress;

	%sketch_cols = 1:floor((size(samples.mean, 2)-1)/S):size(samples.mean, 2);
	sketch_mean = zeros(size(samples.mean, 1), 1);

	S = size(samples.mean, 2);

	for s = 1:S
		[sorted, samples.mean(:, s)] = sort(samples.mean(:, s), 'descend');
		sketch_mean = sketch_mean + sorted / S;
		progress.update(0.7 * s/S);
	end
	for s = 1:S
		samples.mean(samples.mean(:, s), s) = sketch_mean;
		progress.update(0.7 + 0.3 * s/S);
	end
end




fprintf('Summarizing expression data using RMA summarization...\n');

expr = struct;
expr.mean = zeros(length(probesets.probecount), size(samples.mean, 2));

% Reorder the probes so that the probes of each probeset are placed contiguously
% in the vector. The probe index vector therefore has the following form:
% 
%   [ 0 1 2 0 1 0 1 2 3 4 5 6 ]'
% 
% In the example index vector above, we have three probesets with 3, 2 and
% 7 probes, respectively.

probe_reordering = zeros(sum(probesets.probecount), 1);
probe_indices = zeros(sum(probesets.probecount), 1);
k = 1;
for ps = 1:length(probesets.probecount)
	pc = probesets.probecount(ps);
	probe_reordering(k:k+pc-1) = probesets.probes(ps, 1:pc);
	probe_indices(k:k+pc-1) = 0:pc-1;
	k = k + pc;
end

rma_ordered_sample_data = samples.mean(probe_reordering, :);
summary = rmasummary(probe_indices, rma_ordered_sample_data, ...
	'Output', 'natural');
	
% Set the expression of empty probesets to NaN.
expr.mean(:, :) = NaN;
expr.mean(find(probesets.probecount ~= 0), :) = summary;
if isfield(samples, 'meta')
	expr.meta = samples.meta;
	expr.meta.type = probesets.type;
	expr.meta.summarization_method = repmat({'RMA'}, 1, size(samples.mean, 2));
end

