function segments = cgh_segment_cbs(test, ref, probesets, varargin)

global organism;

smooth_window_size = 0;
significance = 0.005;
segment_chr = 1:24;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'segm.*chrom')
		segment_chr = varargin{k+1};
		continue;
	end
	
	if rx(varargin{k}, 'smooth.*win.*size')
		smooth_window_size = varargin{k+1};
		continue;
	end
	
	if rx(varargin{k}, 'significance')
		significance = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

A = test.mean;
B = ref.mean;

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

for chr = segment_chr
	idx = find(probesets.Chromosome == chr);
	ps_in_chr(chr, :) = [min(idx) max(idx)];
end

if smooth_window_size > 0
	for chr = segment_chr
		a = ps_in_chr(chr, 1); b = ps_in_chr(chr, 2);
		logratios(a:b, :) = medfilt2(logratios(a:b, :), [smooth_window_size 1]);
	end
end

segments = struct;
segments.chromosome = cell(length(organism.Chromosomes.Name), size(A, 2));

ps_to_segment = find(ismember(probesets.Chromosome, segment_chr));

for s = 1:size(A, 2)
	% Normalize logratios by moving the highest peak to zero on the x-axis.
	bins = -4:0.05:4;
	n = hist(logratios(:, s), bins);

	bins = bins(2:end-1);
	n = n(2:end-1);

	[~, normal_idx] = max(n);
	normal_level = bins(normal_idx);
	
	logratios(:, s) = logratios(:, s) - normal_level;

	fprintf(1, 'Performing CBS segmentation on sample %s [%d/%d]...\n', ...
		test.meta.sample_id{s}, s, size(A, 2));
	
	% Prepare the data structure that the Matlab CBS function requires.
	cbs_data = struct;
	cbs_data.Chromosome = probesets.Chromosome(ps_to_segment);
	cbs_data.GenomicPosition = probesets.Offset(ps_to_segment);
	cbs_data.Log2Ratio = logratios(ps_to_segment, s);
	
	seg = cghcbs(cbs_data, 'StoppingRule', true, 'Alpha', significance);
	for k = 1:length(seg.SegmentData)
		segments.chromosome{segment_chr(k), s} = seg.SegmentData(k);
	end
end

segments.meta = test.meta;
segments.meta.type = probesets.Type;
segments.meta.segmentation_method = 'Circular binary segmentation';
segments.meta.organism = probesets.Organism;

