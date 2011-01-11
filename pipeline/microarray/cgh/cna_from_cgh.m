function cna = cna_from_cgh(samples, refs, probesets, varargin)

smooth_window_size = 9;
normal_threshold = 0.15;
sample_purity = 1;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'SmoothWindowSize')
		smooth_window_size = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'NormalThreshold')
		normal_threshold = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'SamplePurity')
		sample_purity = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

if mod(smooth_window_size, 2) == 0
	error 'Only odd median window sizes are supported.';
end

if any(size(samples) ~= size(refs))
	error 'The sample and reference matrices must have equal dimensions.';
end

if isfield(samples, 'Mean'), samples = samples.Mean; end
if isfield(refs, 'Mean'), refs = refs.Mean; end

S = size(samples, 2);

cnv = zeros(length(probesets.ProbeCount), S);
for k = 1:length(probesets.ProbeCount)
	probes = probesets.Probes(k, 1:probesets.ProbeCount(k));
	cnv(k, :) = median(samples(probes, :) ./ refs(probes, :), 1);
end

logratios = log2(cnv);

fprintf(1, 'Performing dual window median segmentation...\n');
med_diff = zeros(size(logratios));
for chr = 1:24
	idx = find(probesets.Chromosome == chr);
	a = min(idx); b = max(idx);
	
	med(a:b, :) = medfilt2(logratios(a:b, :), [smooth_window_size 1]);
	
	r = ceil(smooth_window_size / 2);
	d = smooth_window_size;
	med_diff(a+d-1:b-d, :) = med(a+3*r-2:b-r+1, :) - med(a+r-1:b-3*r+2, :);
	med_diff(a:a+d-2, :) = 0;
	med_diff(b-d+1:b, :) = 0;
end

siglevel = 0.0005;

for chr = 1:24
	idx = find(probesets.Chromosome == chr);
	chr_range(chr, :) = [min(idx) max(idx)];
end

for s = 1:S
	mu = mean(med_diff(:, s));
	sigma = std(med_diff(:, s));
	
	p = normcdf(med_diff(:, s), mu, sigma);
	p(p > 0.5) = 1 - p(p > 0.5);
	
	p_bp = (p < siglevel)';
	
	total_bp = 0;
	
	for chr = 1:24
		acmd = abs(med_diff(chr_range(chr, 1):chr_range(chr, 2), s));
		bp = p_bp(chr_range(chr, 1):chr_range(chr, 2));
		
		% Cleanup phase, for contiguous breakpoints we only pick the one with
		% the largest difference of medians.
		run_ends = [find(bp(2:end) ~= bp(1:end-1)) length(bp)];
		run_lengths = diff([0 run_ends]);
		pos = 1;
		for r = 1:length(run_lengths)
			[~, idx] = max(acmd(pos:pos+run_lengths(r)-1));
			bp(pos:pos+run_lengths(r)-1) = false;
			bp(pos+idx-1) = true;
			pos = pos + run_lengths(r);
		end
		
		total_bp = total_bp + sum(bp);
		
		bp = [chr_range(chr, 1) - 1, chr_range(chr, 1) - 1 + find(bp), ...
			chr_range(chr, 2)];   % Breakpoints
		
		for k = 1:length(bp)-1
			seg = bp(k)+1:bp(k+1);
			logratios(seg, s) = median(logratios(seg, s));
		end
	end
	
	
	
	
	% Normalize logratios by moving the highest peak to zero on the x-axis.
	bins = -4:0.05:4;
	n = hist(logratios(:, s), bins);

	bins = bins(2:end-1);
	n = n(2:end-1);

	[~, normal_idx] = max(n);
	normal_level = bins(normal_idx);
	
	logratios(:, s) = logratios(:, s) - normal_level;
	
	
	
	
	% Identify regions of no copy number change.
	run_bounds = [0; find(logratios(2:end, s) ~= logratios(1:end-1, s)); ...
		size(logratios, 1)];
	run_lengths = diff(run_bounds);
	
	seg_ratios = zeros(size(run_lengths));
	
	pos = 1;
	for r = 1:length(run_lengths)
		seg_ratio = logratios(pos, s);
		if abs(seg_ratio) < normal_threshold
			logratios(pos:pos+run_lengths(r)-1, s) = 0;
		end
		pos = pos + run_lengths(r);
	end
end

cna = 2 * 2.^logratios - 2;
cna = cna / sample_purity;
cna(cna < -2) = -2;

% Throw away probes that span too large areas to be informative. This makes
% areas such as centromeres have an undefined copy number, as they should.
for chr = 1:24
	idx = find(probesets.Chromosome == chr);
	offsets = probesets.Offset(idx);
	borders = [offsets(1); round(mean([offsets(1:end-1) offsets(2:end)], ...
		2)); offsets(end)];
	
	too_wide = (borders(2:end) - borders(1:end-1) > 500000);
	cna(idx(too_wide), :) = NaN;
end



