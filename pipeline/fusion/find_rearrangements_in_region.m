
% Author: Matti Annala <matti.annala@tut.fi>

function rearrangements = find_rearrangements_in_region(reads, regions, ...
	min_distance, varargin)

global organism;
chromosomes = organism.Chromosomes;

if ischar(regions), regions = {regions}; end
	
discard_txome_reads = true;
max_mismatches = 1;
anchor_len = 22;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'max.*mismatch')
		max_mismatches = varargin{k+1};
		continue;
	end
	
	if rx(varargin{k}, 'anchor.*len')
		anchor_len = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

S = length(reads.url);

rearrangements = struct;
rearrangements.ReadSequences = {};
rearrangements.Offset = [];




% Parse the genomic regions that the user wishes to calculate rearrangements
% between.
chr = nan(length(regions), 1);
range = nan(length(regions), 2);
region_seq = cell(length(regions), 1);

for k = 1:length(regions)
	tokens = regexpi(regions{k}, '^chr(.+?):\s*(\d+)\s*-\s*(\d+)$', 'tokens');
	if length(tokens) ~= 1, error 'Invalid range specified.'; end
	
	token = tokens{1};
	chr(k) = chromosome_sym2num(token{1});
	range(k, :) = str2double(token(2:3));
	
	if isnan(chr(k)), error 'Invalid chromosome specified.'; end
	
	region_seq{k} = chromosomes.Sequence{chr(k)}(range(k,1):range(k,2));
end
	
	






for s = 1:S
	fprintf(1, 'Analyzing sample %s:\n', reads.meta.sample_id{s});
	
	if ~rx(reads.paired, 'pair')
		fprintf(1, '-> Splitting reads into %dbp start and end anchors...\n',...
			anchor_len);
		split = split_reads(filter(reads, s), anchor_len);
	else
		split = filter(reads, s);
	end
		
	fprintf(1, '-> Aligning paired reads to the region of interest...\n');
	alignments = bowtie_align(split, region_seq, ...
		sprintf('-m1 -v%d', max_mismatches));
	
	fprintf(1, '-> Searching for rearrangement events...\n');
	al = all_alignments(alignments)
	
	read_ids = al.read;
	offsets = al.offset;
	strand = al.strand;
	read_seq = al.sequence;
	target = str2double(al.target);
	
	read_id_idx = containers.Map(read_ids, num2cell(1:length(read_ids)));

	% Last 5p base, first 3p base, 5p strand, 3p strand
	fid = fopen(sprintf('rearrangements_%d.sam', s), 'W');

	for r = 1:length(read_ids)
		id = read_ids{r};
		
		if ~strcmp('/2', id(end-1:end)), continue, end
		
		if read_id_idx.isKey([id(1:end-2) '/1'])
			l = read_id_idx([id(1:end-2) '/1']);
			
			left_offset = offsets(l);
			right_offset = offsets(r);
			
			left_chr = chr(target(l));
			right_chr = chr(target(r));
			
			left_range = range(target(l), :);
			right_range = range(target(r), :);

			if target(l) == target(r) && ...
				abs(right_offset - left_offset) < min_distance
				continue;
			end
			
			left_seq = read_seq{l};
			right_seq = read_seq{r};
			
			flags = 1 + 2 + 64;
			if strand(l) == '-', flags = flags + 16; end
			fprintf(fid, ['r%06d/1\t%d\tchr%s\t%d\t255\t%dM\t=\t' ...
				'%d\t%d\t%s\t*\n'], ...
				r, flags, chromosomes.Name{left_chr}, ...
				left_offset + left_range(1), length(left_seq), ...
				right_offset + right_range(1), ...
				abs(right_offset - left_offset) + length(left_seq), left_seq);
			flags = 1 + 2 + 128;
			if strand(r) == '-', flags = flags + 16; end
			fprintf(fid, ['r%06d/2\t%d\tchr%s\t%d\t255\t%dM\t=\t' ...
				'%d\t%d\t%s\t*\n'], ...
				r, flags, chromosomes.Name{right_chr}, ...
				right_offset + right_range(1), length(right_seq), ...
				left_offset + left_range(1), ...
				abs(right_offset - left_offset) + length(right_seq), right_seq);
		end
	end
	
	fclose(fid);
end


