
% Author: Matti Annala <matti.annala@tut.fi>

function segments = cnvseq_rcbs(cnv, varargin)

significance = 1e-8;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'significance')
		significance = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

S = size(cnv.mean, 2);

tmp = temporary('rcbs');

C = max(cnv.rows.chromosome);
chr = cnv.rows.chromosome;
offset = cnv.rows.position;

segments = struct;
segments.meta = cnv.meta;
segments.chromosome = cell(C, S);

for s = 1:S
	lr = cnv.mean(:, s);
	save(sprintf('%s%d.mat', tmp, s), 'lr', 'chr', 'offset', '-v6');
	
	r_script = sprintf('%s%d.Rscript', tmp, s);

	fid = fopen(r_script, 'W');
	fprintf(fid, [ ...
		'library(R.matlab)\n' ...
		'library(DNAcopy)\n' ...
		sprintf('data = readMat("%s%d.mat")\n', tmp, s) ...
		['cna = CNA(data$lr, data$chr, data$offset, ' ...
			'data.type = "logratio")\n'] ...
		'smoothed.cna = smooth.CNA(cna)\n' ...
		sprintf(['seg.cna = segment(smoothed.cna, alpha=%f, min.width=5, ' ...
			'undo.splits="sdundo", undo.SD=1)\n'], significance) ...
		sprintf(['writeMat("%s%d.mat", chr = seg.cna$output$chrom, ' ...
			'segstart = seg.cna$output$loc.start, ' ...
			'segend = seg.cna$output$loc.end, ' ...
			'segmean = seg.cna$output$seg.mean)\n'], tmp, s)]);
	fclose(fid);
end

fprintf('Performing circular binary segmentation...\n');
unix(['cd ' tmp ' && echo *.mat | parallel -w -n10 ''R CMD BATCH ${x%.mat}.Rscript''']);
	
% FIXME: Maybe grab these from the R analysis instead?
% seg.cna$segRows$start
% seg.cna$segRows$end

for s = 1:S
	results = load(sprintf('%s%d.mat', tmp, s));
	for c = 1:C
		chr_segs = (results.chr == c);
		seg = struct;
		seg.start = results.segstart(chr_segs);
		seg.end = results.segend(chr_segs);
		seg.logratio = results.segmean(chr_segs);
		segments.chromosome{c, s} = seg;
	end
end


