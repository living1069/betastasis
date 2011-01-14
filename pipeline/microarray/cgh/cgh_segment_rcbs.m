
% CGH_SEGMENT_RCBS   Segment aCGH data using R's CBS implementation
%
%    SEGMENTS = CGH_SEGMENT_RCBS(TEST, REF, PROBESETS) calculates copy number 
%    segments using the paired aCGH data in TEST and REF. The mappings from CGH
%    probes to chromosomal loci are provided in the argument PROBESETS. The
%    function uses the CBS implementation of R's DNAcopy package.
%
%    In order to use this function, you must have R installed with the packages
%    DNAcopy, R.matlab and any package dependencies of these.

% Author: Matti Annala <matti.annala@tut.fi>

function segments = cgh_segment_rcbs(test, ref, probesets, varargin)

global organism;

normal_threshold = 0.2;
sample_purity = 0.7;
significance = 1e-8;
detect_gender = false;

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
	
	%if strcmpi(varargin{k}, 'DetectGender')
	%	detect_gender = varargin{k+1};
	%	continue;
	%end
	
	error('Unrecognized option "%s".', varargin{k});
end

logratios = cgh_to_logratios(test, ref, probesets, 'Smooth', 0);
chr = probesets.Chromosome;
offset = probesets.Offset;

S = size(logratios, 2);

tmp = ptemp

r_script = [tmp '.Rscript'];

fid = fopen(r_script, 'W');
fprintf(fid, [ ...
	'library(R.matlab)\n' ...
	'library(DNAcopy)\n' ...
	'data = readMat("' [tmp '.mat'] '")\n' ...
	'cna = CNA(data$lr, data$chr, data$offset, data.type = "logratio")\n' ...
	'smoothed.cna = smooth.CNA(cna)\n' ...
	'segmented.cna = segment(smoothed.cna, verbose = 1)\n' ...
	'writeMat("' [tmp '.mat'] '", chr = segmented.cna$output$chrom, segstart = segmented.cna$output$loc.start, segend = segmented.cna$output$loc.end, segmean = segmented.cna$output$seg.mean)\n']);
fclose(fid);


ploidy = cgh_ploidy(test, ref);


segments = struct;
segments.Chromosome = cell(length(organism.Chromosomes.Name), S);

progress = Progress;

for s = 1:S
	lr = logratios(:, s);
	save([tmp '.mat'], 'lr', 'chr', 'offset', '-v6');
	
	[status, out] = unix(sprintf('/worktmp/R/bin/R CMD BATCH %s ~/out.txt', ...
		r_script));
	if status ~= 0, error 'R CBS segmentation failed.'; end
	
	results = load([tmp '.mat']);
	
	for c = 1:length(organism.Chromosomes.Name)
		chr_segs = (results.chr == c);
		seg = struct;
		
		if ~isnan(ploidy(c, s))
			seg.Start = results.segstart(chr_segs);
			seg.End = results.segend(chr_segs);
		
			lr = results.segmean(chr_segs);
			seg.CNA = ploidy(c, s) * 2.^lr - ploidy(c, s);
			seg.CNA(seg.CNA < - 2) = -2;
			
			% Filter out segments below the normal threshold.
			seg.CNA(abs(seg.CNA) < normal_threshold) = 0;
		end
		
		segments.Chromosome{c, s} = seg;
	end
		
	progress.update(s / S);
end

segments.Meta = test.Meta;
segments.Meta.Ref = ref.Meta;
segments.Meta.Organism = probesets.Organism;
segments.Meta.Type = 'Copy number segments';
segments.Meta.SegmentationMethod = repmat({'CBS (R)'}, S, 1);
segments.Meta.SamplePurity = sample_purity;
segments.Meta.NormalThreshold = normal_threshold;
segments.Meta.Ref = rmfield(segments.Meta.Ref, 'Type');
segments.Meta.Ref = rmfield(segments.Meta.Ref, 'Platform');

