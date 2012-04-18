
% Author: Matti Annala <matti.annala@tut.fi>

function stats = plot_seq_reads(reads, sequence, varargin)
	
global organism;
exons = organism.Exons;

fragments = [];

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Fragments')
		fragments = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

sample_ids = reads.Meta.Sample.ID;
S = length(sample_ids);

stats = struct;
stats.Samples = sample_ids;
stats.Reads = cell(1, S);

% We iterate through the samples first, and then we iterate through their
% technical replicates.
for s = 1:S
	strands = '';
	offsets = [];
	sequences = {};
	
	index.Sequence = { sequence };
	
	al = align_reads(filter_query(reads, s), index, ...
		'MaxMismatches', 2, 'Columns', 'strand,offset,sequence', varargin{:});
	
	stats.TotalReads(s) = sum(al.TotalReads);
	stats.AlignedReads(s) = sum(al.AlignedReads);

	strands = al.Strand;
	offsets = al.Offset;
	sequences = al.Sequence;
	
	
	
	
	% Calculate arrow positions.
	x = zeros(7, length(offsets));
	y = zeros(7, length(offsets));
	
	offsets = offsets(randperm(length(offsets)));
	
	read_coverage = false(1000, length(sequence));
	
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
		dx = length(sequence) / 800;
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
	figure; hold all; xlim([0 length(sequence)]); ylim([-9 ymax]);
	patch(x, y, 'k');
	
	if ~isempty(fragments)
		for m = 1:size(tx_exons, 1)
			left = exon_pos(m, 1);
			right = exon_pos(m, 2);
			fill([left right right left], [-5 -5 -2 -2], [.9 .9 .9]);
		end
	else
		fill([0 length(sequence) length(sequence) 0], ...
			[-5 -5 -2 -2], [.9 .9 .9]);
	end

	xlabel('Read offset in transcript sequence');
	ylabel('Number of overlapping reads');
	saveas(gcf, sprintf('reads_%s.pdf', sample_ids{s}));
end

