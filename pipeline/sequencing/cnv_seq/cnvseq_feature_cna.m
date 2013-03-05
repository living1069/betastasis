
% Author: Matti Annala <matti.annala@tut.fi>

function fcna = cnvseq_feature_cna(cna, features, varargin)

min_region = 0;

% FIXME: Remove this once the naming is unified.
if ~isfield(features, 'chromosome') && isfield(features, 'Chromosome')
	features.chromosome = features.Chromosome;
	features.position = features.Position;
end

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'min.*region')
		min_region = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

S = size(cna.mean, 2);
F = length(features.chromosome);

fcna = struct;
fcna.meta = cna.meta;
fcna.rows = features;
fcna.mean = nan(F, S);

% If chromosomes are expressed as a cell array, convert them to numbers
if iscellstr(features.chromosome)
	features.chromosome = chromosome_sym2num(features.chromosome);
end

% Make sure that all feature positions are expressed lower coordinate first
fpos = features.position;
wrong = fpos(:, 2) < fpos(:, 1);
fpos(wrong, :) = [fpos(wrong, 2), fpos(wrong, 1)];

% Extend features to a minimum size of "min_region"
too_small = fpos(:, 2) - fpos(:, 1) < min_region;
d = min_region - (fpos(:, 2) - fpos(:, 1));
fpos(too_small, 1) = fpos(too_small, 1) - round(d(too_small) / 2);
fpos(too_small, 2) = fpos(too_small, 2) + round(d(too_small) / 2);

for f = find(all(~isnan(fpos), 2))'
	% Find all windows that target the feature we're interested in.
	idx = find(cna.rows.chromosome == features.chromosome(f) & ...
		cna.rows.position >= fpos(f, 1) & cna.rows.position <= fpos(f, 2));
	if isempty(idx), continue, end
		
	fcna.mean(f, :) = nanmedian(cna.mean(idx, :), 1);
end

