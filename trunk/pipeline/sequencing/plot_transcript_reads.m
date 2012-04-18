
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

function stats = plot_transcript_reads(reads, transcript, image_prefix, ...
	varargin)
	
global organism;
exons = organism.Exons;

if ischar(transcript)
	tx_idx = transcript_idx(transcript);
elseif isnumeric(transcript)
	tx_idx = transcript;
end

tx_exons = organism.Transcripts.Exons{tx_idx};
tx_exon_pos = organism.Transcripts.ExonPos{tx_idx};

sample_id = reads.Meta.Sample.ID;

if isfield(sample_id, 'ID')
	[~, uniq_samples] = unique(sample_id);
else
	uniq_samples = 1:length(reads.Raw);
end

S = length(uniq_samples);

tx_seq = organism.Transcripts.Sequence{tx_idx};

stats.Samples = sample_id(uniq_samples);

stats.TotalReads = zeros(1, S);
stats.AlignedReads = zeros(1, S);

stats.ExonReads = zeros(length(tx_exons), S);
stats.JunctionReads = zeros(length(tx_exons) - 1, S);

% We iterate through the samples first, and then we iterate through their
% technical replicates.
for s = 1:S
	fprintf(1, 'Plotting reads for sample %s...\n', sample_id{s});
	replicates = find(strcmp(sample_id{s}, sample_id));
	
	strands = '';
	offsets = [];
	sequences = {};
	
	if 1
		index.Name = { organism.Transcripts.Name{tx_idx} };
		index.Sequence = { organism.Transcripts.Sequence{tx_idx} };
	else
		index = 'transcripts';
	end
	
	% If we have technical replicates (multiple sequenced channels) from a
	% particular sample, we combine the alignments from all the replicates and
	% then plot them together for that sample.
	for r = replicates'
		if length(replicates) > 1
			fprintf(1, '-> Channel %s...\n', reads.Meta.Sample.Filename{r});
		end
		
		al = align_reads(filter_query(reads, r), index, ...
			'MaxMismatches', 2, ...
			'Columns', 'target,strand,offset,sequence', varargin{:});
			
		stats.TotalReads(s) = stats.TotalReads(s) + sum(al.TotalReads);
		stats.AlignedReads(s) = stats.AlignedReads(s) + sum(al.AlignedReads);

		keep = strcmp(organism.Transcripts.Name{tx_idx}, al.Target);
		
		strands = cat(1, strands, al.Strand(keep));
		offsets = cat(1, offsets, al.Offset(keep));
		sequences = cat(1, sequences, al.Sequence(keep));
	end
	
	
	
	
	% Calculate the number of reads aligned to exons and exon junctions.
	[ordered_offsets, order] = sort(offsets);
	ex = 1;
	ex_pos = tx_exon_pos(ex, :);
	for m = 1:length(ordered_offsets)
		while 1
			read_pos = [ordered_offsets(m), ...
				ordered_offsets(m) + length(sequences{order(m)}) - 1];
			if read_pos(1) > ex_pos(2)
				ex = ex + 1;
				ex_pos = tx_exon_pos(ex, :);
			elseif read_pos(2) > ex_pos(2)
				stats.JunctionReads(ex, s) = stats.JunctionReads(ex, s) + 1;
				break;
			else
				stats.ExonReads(ex, s) = stats.ExonReads(ex, s) + 1;
				break;
			end
		end
	end
	
	
	
	
	
	
	
	% Calculate arrow positions.
	x = zeros(7, length(offsets));
	y = zeros(7, length(offsets));
	
	offsets = offsets(randperm(length(offsets)));
	
	read_coverage = false(1000, length(tx_seq));
	
	highest_level = -Inf;

	for m = 1:length(offsets)
		span = offsets(m):offsets(m)+length(sequences{m})-1;
		
		c = 1;
		while any(read_coverage(c, span))
			c = c + 1;
		end
		
		read_coverage(c, span) = true;
		
		level = c;
		highest_level = max(highest_level, level);
		
		len = length(sequences{m});
		dx = length(tx_seq) / 800;
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
	


	% Render the quiver plot.
	ymax = max(150, highest_level + 10);
	figure; hold all; xlim([0 length(tx_seq)]); ylim([-9 ymax]);
	patch(x, y, 'k');
	
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
		
	title(sprintf('Quiver plot of reads for transcript %s%s\nin sample %s', ...
		organism.Transcripts.Name{tx_idx}, gene_name, sample_id{s}), ...
		'Interpreter', 'none');
	xlabel('Read offset in transcript sequence');
	ylabel('Number of overlapping reads');
	saveas(gcf, [image_prefix '_' num2str(s) '.pdf']);
end

