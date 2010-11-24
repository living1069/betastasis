function segments = cgh_segment_cbs(samples, refs, probesets, varargin)

global organism;

smooth_window_size = 9;

normal_threshold = 0.2;
sample_purity = 0.7;
significance = 0.005;
drop_sex_chromosomes = true;
show_segments = false;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'SmoothWindowSize')
		smooth_window_size = varargin{k+1};
		continue;
	end
	
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
	
	if strcmpi(varargin{k}, 'SexChromosomes')
		drop_sex_chromosomes = ~varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

A = samples.Mean;
B = refs.Mean;

if any(size(A) ~= size(B))
	error 'The sample and reference matrices must have equal dimensions.';
end

S = size(A, 2);

cnv = zeros(length(probesets.ProbeCount), S);
for k = 1:length(probesets.ProbeCount)
	probes = probesets.Probes(k, 1:probesets.ProbeCount(k));
	cnv(k, :) = median(A(probes, :), 1) ./ median(B(probes, :), 1);
end

logratios = log2(cnv);

for chr = 1:24
	idx = find(probesets.Chromosome == chr);
	chr_range(chr, :) = [min(idx) max(idx)];
end

for chr = 1:24
	a = chr_range(chr, 1); b = chr_range(chr, 2);
	logratios(a:b, :) = medfilt2(logratios(a:b, :), [smooth_window_size 1]);
end

for s = 1:size(A, 2)
	% Normalize logratios by moving the highest peak to zero on the x-axis.
	bins = -4:0.05:4;
	n = hist(logratios(:, s), bins);

	bins = bins(2:end-1);
	n = n(2:end-1);

	[~, normal_idx] = max(n);
	normal_level = bins(normal_idx);
	
	logratios(:, s) = logratios(:, s) - normal_level;
end

logratio_tmp = ptemp
fid = fopen(logratio_tmp, 'W');
for k = 1:length(probesets.ProbeCount)
	fprintf(fid, 'chr%s\t%d', organism.Chromosomes.Name{ ...
		probesets.Chromosome(k)}, probesets.Offset(k));
	fprintf(fid, '\t%f', logratios(k, :));
	fprintf(fid, '\n');
end
fclose(fid);
return;


segments = struct;
segments.Chromosome = cell(length(organism.Chromosomes.Name), size(A, 2));

for s = 1:size(A, 2)
	fprintf(1, 'Performing CBS segmentation on sample %s [%d/%d]...\n', ...
		samples.Meta.Sample.ID{s}, s, size(A, 2));
	
	N = length(probesets.ProbeCount);
		
	% Prepare the data structure that the Matlab CBS function requires.
	cbs_data = struct;
	cbs_data.Chromosome = probesets.Chromosome;
	cbs_data.GenomicPosition = probesets.Offset;
	cbs_data.Log2Ratio = logratios;
	
	seg = cghcbs(cbs_data, 'StoppingRule', true, 'Alpha', significance);
	for k = 1:length(seg.SegmentData)
		segments.Chromosome{k, s} = seg.SegmentData(k);
	end
end

segments.Meta = samples.Meta;
segments.Meta.Type = probesets.Type;
segments.Meta.SegmentationMethod = 'Circular binary segmentation';
segments.Meta.Organism = probesets.Organism;
segments.Meta.Platform = repmat({segments.Meta.Platform}, ...
	size(samples.Mean, 2), 1);

	
	
	
	
	
	





