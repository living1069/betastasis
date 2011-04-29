
% PLOT_TRANSCRIPT_READS   Draw a quiver plot of alignments for an RNA transcript
% 
%    PLOT_TRANSCRIPT_READS(READS, TX, IMAGE_PREFIX, ...) aligns the
%    sequencing samples READS against the sequence of the RNA transcript
%    TX, where TX is a transcript name (e.g. 'NM_005560') or index into 
%    organism.Transcripts. IMAGE_PREFIX is used for choosing a filename
%    for the PDF alignment plots that the function produces. The full filename
%    for the PDF files will be of the form IMAGE_PREFIX_1.pdf, where the number
%    at the end is the sample index in READS.
%    
%    The alignment step can be customized by providing additional options
%    arguments; see HELP ALIGN_READS for a full list of available options.

% Author: Matti Annala <matti.annala@tut.fi>

function [] = plot_transcript_reads(reads, transcript, image_prefix, ...
	varargin)
	
global organism;

if ischar(transcript)
	tx_idx = transcript_idx(transcript);
elseif isnumeric(transcript)
	tx_idx = transcript;
end

S = length(reads.Raw);

for k = 1:S
	tx_seq = organism.Transcripts.Sequence{tx_idx};
	
	al = align_reads(filter_query(reads, k), {tx_seq}, varargin{:}, ...
		'MaxMismatches', 2, 'Columns', 'strand,offset,sequence');

	strands = al.Strand;
	offsets = al.Offset;
	sequences = al.Sequence;
	
	read_count_dist = zeros(1, length(tx_seq));
	x = zeros(7, length(offsets));
	y = zeros(7, length(offsets));

	for m = 1:length(offsets)
		span = offsets(m):offsets(m)+length(sequences{m})-1;
		read_count_dist(span) = read_count_dist(span) + 1;
		
		level = max(read_count_dist(span));
		len = length(sequences{m});
		dx = length(read_count_dist) / 800;
		dy = 150 / 800; 
		
		arrow_x = [0; 0; len-dx; len-dx; len; len-dx; len-dx];
		arrow_y = [-dy/2; dy/2; dy/2; dy; 0; -dy; -dy/2];
		
		if strcmp(strands(m), '+')
			x(:, m) = arrow_x + double(offsets(m));
		else
			x(:, m) = (len - arrow_x) + double(offsets(m));
		end
		y(:, m) = arrow_y + level;
    end
	
%	figure; bar(1:length(read_count_dist), read_count_dist);
%	title(sprintf('Read distribution in sample %s', reads.Meta.Sample.ID{k}),...
%		'Interpreter', 'none');
%	xlabel('Feature sequence offset'); ylabel('Number of overlapping reads');
%	saveas(gcf, [image_prefix '_' num2str(k) '.pdf']);
	
%	figure; bar(1:length(read_count_dist), ...
%		conv(read_count_dist, ones(1, 500) / 500, 'same'));
%	title(sprintf('Smoothed read distribution in sample %s', ...
%		reads.Meta.Sample.ID{k}));
%	xlabel('Feature sequence offset'); ylabel('Number of overlapping reads');
%	saveas(gcf, [image_prefix '_smoothed_' num2str(k) '.pdf']);
	
	figure; hold all; xlim([0 length(tx_seq)]); ylim([-9 150]);
	patch(x, y, 'k');
	
	tx_exons = organism.Transcripts.Exons{tx_idx};
	exon_pos = organism.Transcripts.ExonPos{tx_idx};
	for m = 1:length(tx_exons)
		left = exon_pos(m, 1);
		right = exon_pos(m, 2);
		fill([left right right left], [-5 -5 -2 -2], [.9 .9 .9]);
	end
	
	if isempty(tx_exons)
		fill([0 length(tx_seq) length(tx_seq) 0], [-5 -5 -2 -2], [.9 .9 .9]);
	end
	
	cds = organism.Transcripts.CDS(tx_idx, :);
	if ~any(isnan(cds))
		fill([cds(1) cds(2) cds(2) cds(1)], [-4 -4 -3 -3], [.5 .5 .5]);
	end
	
	gene_name = [' (' organism.Genes.Name{ ...
		organism.Transcripts.Gene(tx_idx)} ')'];
		
	if isfield(reads.Meta.Sample, 'ID')
		sample_id = reads.Meta.Sample.ID{k};
	else
		sample_id = reads.Meta.Sample.Filename{k};
	end

	title(sprintf('Quiver plot of reads for transcript %s%s\nin sample %s', ...
		organism.Transcripts.Name{tx_idx}, gene_name, sample_id), ...
		'Interpreter', 'none');
	xlabel('Read offset in transcript sequence');
	ylabel('Number of overlapping reads');
	saveas(gcf, [image_prefix '_quiver_' num2str(k) '.pdf']);
end

