
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

significance = 0.01;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'significance')
		significance = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

logratios = cgh_to_logratios(test, ref, probesets, 'Smooth', 0);

S = size(logratios, 2);

tmp = temporary('rcbs');
r_script = [tmp '.Rscript'];

fid = fopen(r_script, 'W');
fprintf(fid, [ ...
	'library(R.matlab)\n' ...
	'library(DNAcopy)\n' ...
	'data = readMat("' [tmp '.mat'] '")\n' ...
	'cna = CNA(data$lr, data$chr, data$offset, data.type = "logratio")\n' ...
	'smoothed.cna = smooth.CNA(cna)\n' ...
	sprintf('seg.cna = segment(smoothed.cna, alpha=%f)\n', significance), ...
	'writeMat("' [tmp '.mat'] '", chr = seg.cna$output$chrom, segstart = seg.cna$output$loc.start, segend = seg.cna$output$loc.end, segmean = seg.cna$output$seg.mean)\n']);
fclose(fid);

segments = struct;
segments.chromosome = cell(length(organism.Chromosomes.Name), S);

progress = Progress;

for s = 1:S
	fprintf('Performing circular binary segmentation for sample %s...\n', ...
		test.meta.sample_id{s});
		
	chr = probesets.Chromosome;
	offset = probesets.Offset;
	lr = logratios(:, s);
	save([tmp '.mat'], 'lr', 'chr', 'offset', '-v6');
	
	[status, out] = unix(sprintf('R CMD BATCH %s ~/out.txt', r_script));
	if status ~= 0, error 'R CBS segmentation failed.'; end
	
	results = load([tmp '.mat']);
	
	for c = 1:length(organism.Chromosomes.Name)
		chr_segs = (results.chr == c);
		seg = struct;
		seg.start = results.segstart(chr_segs);
		seg.end = results.segend(chr_segs);
		seg.logratio = results.segmean(chr_segs);
		segments.chromosome{c, s} = seg;
	end
		
	progress.update(s / S);
end

segments.meta = test.meta;
segments.meta.type = 'aCGH segments (logratio)';
segments.meta.segmentation_method = repmat({'CBS (R)'}, S, 1);

