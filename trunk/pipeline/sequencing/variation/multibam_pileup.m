
function [] = multibam_pileup(locus, varargin)

global organism;
chromosomes = organism.Chromosomes;

max_read_len = 200;
bam_root = '.';
highlight = [];      % Bases to highlight in output
pileup_dir = '.';

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'bam.*root')
		bam_root = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'highlight')
		highlight = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'pileup.*dir')
		pileup_dir = varargin{k+1}; continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

files = find_files('.bam$', bam_root);


% Parse the genomic coordinates for which we want to render the pileup.
tokens = regexp(locus, '(chr.+?):(\d+)-(\d+)', 'tokens');
if ~isempty(tokens)
	token = tokens{1}; chr = token{1}; offset = str2double(token{2});
	region = [str2double(token{2}) str2double(token{3})];
	pileup_file = sprintf('%s:%d-%d_pileup.txt', chr, region(1), region(2));
end

tokens = regexp(locus, '(chr.+?):(\d+)$', 'tokens');
if ~isempty(tokens)
	token = tokens{1}; chr = token{1}; offset = str2double(token{2});
	radius = 50;
	region = [offset-radius offset+radius];
	if isempty(highlight), highlight = offset; end
	pileup_file = sprintf('%s:%d_pileup.txt', chr, offset);
end





% This is the padded region containing extra flanks that guarantee that we get
% also all reads that only partially overlap with the region of interest.
pad_region = [region(1)-max_read_len region(2)+max_read_len];
highlight = highlight - pad_region(1) + 1;

ref_seq = chromosomes.Sequence{chromosome_sym2num(chr)} ...
	(pad_region(1):pad_region(2));
	
if ~strcmp(pileup_dir, '.'), [~, ~] = mkdir(pileup_dir); end
fid = fopen([pileup_dir '/' pileup_file], 'W');

for f = 1:length(files)
	fprintf(fid, '%s:\n%s\n', files{f}, repmat('=', 1, length(files{f})+1));
	
	al = bamread(files{f}, chr, pad_region, 'Full', true);
	%fprintf(fid, '%d\t%d\t%s\t%s\n', al(1).Position, al(1).Flag, ...
	%	al(1).CigarString, al(1).Sequence);
	rel_pos = [al.Position] - pad_region(1) + 1;
	
	% The cigar2align() function fails utterly if option "GapsInRef" is
	% switched on. So we turn it off and handle things manually instead.
	mal = cigar2align({al.Sequence, ref_seq}, ...
		{al.CigarString, sprintf('%dM', length(ref_seq))}, ...
		'Start', [rel_pos, 1]);
	
	% Only highlighted bases should be shown in uppercase.
	non_hl = setdiff(1:(pad_region(2)-pad_region(1)+1), highlight);
	mal(:, non_hl) = lower(mal(:, non_hl));
	mal(:, highlight) = upper(mal(:, highlight));
	mal(end, :) = upper(mal(end, :));    % Refseq in uppercase
		
	S = size(mal, 1);     % Number of aligned reads + reference sequence.
		
	inserts = repmat({''}, size(mal));
	for k = 1:length(al)
		cigar = al(k).CigarString;
		
		% As we go along the CIGAR string, we track the position in the matrix
		% and along the read sequence.
		ref_offset = rel_pos(k);
		read_offset = 1;
		
		tokens = regexp(cigar, '(\d+)(.)', 'tokens');
		for t = 1:length(tokens)
			clen = str2double(tokens{t}{1}); ctype = tokens{t}{2};
			if strcmp(ctype, 'D')
				ref_offset = ref_offset + clen;
			elseif strcmp(ctype, 'I')
				inserts{k, ref_offset} = ...
					al(k).Sequence(read_offset:read_offset+clen-1);
				read_offset = read_offset + clen;
			else
				ref_offset = ref_offset + clen;
				read_offset = read_offset + clen;
			end
		end
	end
	
	total_gaps = 0;
	imal = repmat(' ', S, 0);
	for k = 1:size(inserts, 2)
		max_ilen = max(cellfun(@length, inserts(:, k)));
		if max_ilen > 0
			gap = repmat('-', S, max_ilen);
			
			% Clear the '-' on lines where reads have already ended.
			gap(imal(:, end) == ' ', :) = ' ';
			
			for s = 1:S
				if isempty(inserts{s,k}), continue, end
				gap(s, 1:length(inserts{s,k})) = inserts{s,k};
			end
			imal = [imal gap];
			total_gaps = total_gaps + size(gap, 2);
		end
		imal = [imal mal(:, k)];
	end
	
	% Take into account that the refseq length is increased by gaps.
	%line_lengths = cellfun(@length, cellstr(imal));
	%imal = imal(line_lengths <= length(ref_seq) + total_gaps, :);
	
	% Throw away empty lines and trim to the desired region.
	imal = imal(:, max_read_len+1:end-max_read_len);
	imal = imal(any(imal ~= ' ', 2), :);
	imal = cellstr(imal);
	
	fprintf(fid, '%s\n', imal{:});
	fprintf(fid, '\n');
end

fclose(fid);

