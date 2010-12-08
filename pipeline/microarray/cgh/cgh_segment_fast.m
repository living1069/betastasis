
% CGH_SEGMENT_FAST   Segment array CGH data using multiscale difference of means
%
%    SEGS = CGH_SEGMENT_FAST(TEST, REF, PROBESETS) calculates copy number 
%    segments using the paired aCGH data in TEST and REF. The mappings from CGH
%    probes to chromosomal loci are provided in the argument PROBESETS. The
%    algorithm performs multiscale mean filtering on the CGH ratios, and then
%    detects segment breakpoints by looking for edges in the filtered signal.
% 
%    CGH_SEGMENT_FAST(..., 'Significance', ALPHA) specifies the significance
%    level ALPHA at which breakpoints will be called. A high significance 
%    threshold gives high sensitivity, a low threshold high specificity.
%    Default significance threshold is 1e-8.
%
%    CGH_SEGMENT_FAST(..., 'SamplePurity', PURITY) specifies the average sample
%    purity. All calculated copy number alterations will be scaled up by
%    1/PURITY. Default sample purity is 0.7.
%
%    CGH_SEGMENT_FAST(..., 'NormalThreshold', NT) specifies a CNA threshold
%    for calling unaltered chromosomal segments. That is, if for any segment
%    |CNA| < NT, the segment will be marked as unaltered. Default is 0.2.

% Author: Matti Annala <matti.annala@tut.fi>

function segments = cgh_segment_fast(test, ref, probesets, varargin)

global organism;

normal_threshold = 0.2;
sample_purity = 0.7;
significance = 1e-8;
detect_gender = false;
smooth_window_size = 7;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'NormalThreshold')
		normal_threshold = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Significance')
		significance = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'SamplePurity')
		sample_purity = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'SmoothWindowSize')
		smooth_window_size = varargin{k+1};
		continue;
	end
	
	%if strcmpi(varargin{k}, 'DetectGender')
	%	detect_gender = varargin{k+1};
	%	continue;
	%end
	
	error('Unrecognized option "%s".', varargin{k});
end

% We apply a slight median filtering to the data before the multiscale
% analysis that is based on mean filtering. The median filtering helps remove
% artifacts.
logratios = cgh_to_logratios(test, ref, probesets, ...
	'Smooth', smooth_window_size);

S = size(logratios, 2);

chr_range = nan(length(organism.Chromosomes.Name), 2);
for chr = 1:length(organism.Chromosomes.Name)
	idx = find(probesets.Chromosome == chr);
	if ~isempty(idx)
		chr_range(chr, :) = [min(idx) max(idx)];
	end
end



fprintf(1, 'Checking sample ploidy information...\n');

genders_match = false(S, 1);
if isfield(test.Meta, 'Patient') && isfield(test.Meta.Patient, 'Gender') && ...
	isfield(ref.Meta, 'Patient') && isfield(ref.Meta.Patient, 'Gender')
	
	genders_known = ~strcmpi('-', test.Meta.Patient.Gender) & ...
		~strcmpi('-', ref.Meta.Patient.Gender);
	genders_match = strcmpi(test.Meta.Patient.Gender, ...
		ref.Meta.Patient.Gender) | genders_known;
end

if ~all(genders_match)
	fprintf(1, ['WARNING: Skipping sex chromosomes in %d samples with ' ...
	            'mismatched or missing gender information.\n'], ...
				sum(~genders_match));
end

ploidy = nan(length(organism.Chromosomes.Name), S);
for s = 1:S
	if ~genders_match(s) || strcmp(test.Meta.Patient.Gender{s}, '-')
		ploidy(:, s) = [2*ones(1, 22) NaN NaN NaN]';
	elseif strcmpi(test.Meta.Patient.Gender{s}, 'Male')
		ploidy(:, s) = [2*ones(1, 22) 1 1 NaN]';
	elseif strcmpi(test.Meta.Patient.Gender{s}, 'Female')
		ploidy(:, s) = [2*ones(1, 22) 2 0 NaN]';
	else
		ploidy(:, s) = [2*ones(1, 22) NaN NaN NaN]';
		fprintf(1, 'Ploidy of gender "%s" is unknown.\n', ...
			test.Meta.Patient.Gender{s});
	end
end

cna = zeros(size(logratios));
for s = 1:S
	for chr = 1:length(organism.Chromosomes.Name)
		cr = chr_range(chr, 1):chr_range(chr, 2);
		if any(isnan(cr)), continue, end
		
		% Here we mask regions of unknown ploidy with zeros, as using NaNs
		% would make it more difficult to implement the segmentation algorithm.
		% We replace the CNA for these regions with NaNs later in the function.
		if isnan(ploidy(chr, s))
			cna(cr, s) = 0;
		else
			cna(cr, s) = ploidy(chr, s) * 2.^logratios(cr, s) - ploidy(chr, s);
		end
	end
end

cna = cna / sample_purity;
cna(cna < -2) = -2;




fprintf(1, 'Detecting breakpoints with significance threshold %.1x...\n', ...
	significance);

all_bp = false(size(cna));

% Generate scales that are all odd.
scales = 7 * 1.5.^(0:13);
scales = round(scales);
scales = scales - (1 - mod(scales, 2));

progress = Progress;

% Run a multiscale breakpoint analysis.
for k = 1:length(scales)
	scale = scales(k);
	
	med_diff = zeros(size(cna));
	for chr = 1:length(organism.Chromosomes.Name)
		a = chr_range(chr, 1); b = chr_range(chr, 2);
		if isnan(a) || isnan(b), continue, end
		
		med(a:b, :) = conv2(cna(a:b, :), ones(scale, 1) / scale, 'same');
		%med(a:b, :) = medfilt2(cna(a:b, :), [scale 1]);
		
		r = ceil(scale / 2); d = scale;
		range = a+d-1:b-d;
		med_diff(range, :) = med(range + r, :) - med(range - (r - 1), :);
		%med_diff(a:a+d-2, :) = 0;
		%med_diff(b-d+1:b, :) = 0;
	end

	for s = 1:S
		%mu = mean(med_diff(:, s))
		mu = 0;
		
		up = (med_diff(:, s) >= mu);
		down = (med_diff(:, s) < mu);
		
		sigma_up = sqrt(sum((med_diff(up, s) - mu).^2) / (sum(up) - 1));
		sigma_down = sqrt(sum((med_diff(up, s) - mu).^2) / (sum(down) - 1));
		
		p = zeros(size(med_diff, 1), 1);
		p(up) = 1 - normcdf(med_diff(up, s), mu, sigma_up);
		p(down) = normcdf(med_diff(down, s), mu, sigma_down);
		
		bp = (p < significance);
		
		for chr = 1:length(organism.Chromosomes.Name)
			cr = chr_range(chr, 1):chr_range(chr, 2);
			if any(isnan(cr)), continue, end

			acmd = abs(med_diff(cr, s));
			cbp = bp(cr);
			
			% Cleanup phase, for contiguous breakpoints we only pick the one
			% with the largest difference of medians.
			run_ends = [find(cbp(2:end) ~= cbp(1:end-1)); length(cbp)];
			run_lengths = diff([0; run_ends]);
			pos = 1;
			for r = 1:length(run_lengths)
				[~, idx] = max(acmd(pos:pos+run_lengths(r)-1));
				cbp(pos:pos+run_lengths(r)-1) = false;
				cbp(pos+idx-1) = true;
				pos = pos + run_lengths(r);
			end
			
			all_bp(cr, s) = all_bp(cr, s) | cbp;
		end
	end
	
	progress.update(k / length(scales));
end
		
for s = 1:S
	for chr = 1:length(organism.Chromosomes.Name)
		a = chr_range(chr, 1); b = chr_range(chr, 2);
		if isnan(a) || isnan(b), continue, end

		all_cbp = all_bp(a:b, s);
		cbp = [a - 1; a - 1 + find(all_cbp); b];   % Breakpoints
		
		for k = 1:length(cbp)-1
			seg = cbp(k)+1:cbp(k+1);
			cna(seg, s) = median(cna(seg, s));
		end
	end
end

% Identify regions of no copy number change.
cna(abs(cna) < normal_threshold) = 0;

% Mask out genomic regions of unknown ploidy.
for s = 1:S
	for chr = 1:length(organism.Chromosomes.Name)
		if isnan(ploidy(chr, s))
			a = chr_range(chr, 1); b = chr_range(chr, 2);
			if isnan(a) || isnan(b), continue, end

			cna(a:b, s) = NaN;
		end
	end
end




% Throw away probes that span too large areas to be informative. This makes
% areas such as centromeres have an undefined copy number, as they should.
for chr = 1:length(organism.Chromosomes.Name)
	cr = chr_range(chr, 1):chr_range(chr, 2);
	if any(isnan(cr)), continue, end

	offsets = probesets.Offset(cr);
	borders = [offsets(1); round(mean([offsets(1:end-1) offsets(2:end)], ...
		2)); offsets(end)];
	
	too_wide = (borders(2:end) - borders(1:end-1) > 500000);
	cna(cr(too_wide), :) = NaN;
end




% Build the segments.
segments = struct;
segments.Meta = test.Meta;
segments.Meta.Ref = ref.Meta;
segments.Meta.Organism = probesets.Organism;
segments.Meta.Type = 'Copy number segments';
segments.Meta.SegmentationMethod = ...
	repmat({'Multiscale difference of medians'}, S, 1);
segments.Meta.Ref = rmfield(segments.Meta.Ref, 'Type');
segments.Meta.Ref = rmfield(segments.Meta.Ref, 'Platform');

segments.Chromosome = cell(length(organism.Chromosomes.Name), S);

for s = 1:S
	for chr = 1:length(organism.Chromosomes.Name)
		seg = struct;
		seg.Start = [];
		seg.End = [];
		seg.CNA = [];

		cr = chr_range(chr, 1):chr_range(chr, 2);
		if any(isnan(cr))
			segments.Chromosome{chr, s} = seg;
			continue;
		end

		offsets = probesets.Offset(cr);
		ccna = cna(cr, s);
		
		N = length(cr);
		run_bounds = [0; find((ccna(2:end) ~= ccna(1:end-1)) & ...
			~(isnan(ccna(2:end)) & isnan(ccna(1:end-1)))); N];
		run_cna = ccna(run_bounds(1:end-1)+1);
		
		borders = round(mean([offsets(1:N-1) offsets(2:N)], 2));
		borders = [offsets(1); borders; offsets(N)];
		
		run_borders = borders(run_bounds + 1);
		
		for p = 1:length(run_cna)
			seg.Start(end+1) = run_borders(p) + 1;
			seg.End(end+1) = run_borders(p+1);
			seg.CNA(end+1) = run_cna(p);
		end
		
		segments.Chromosome{chr, s} = seg;
	end
end

