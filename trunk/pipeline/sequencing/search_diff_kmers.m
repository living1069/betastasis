
function expr = search_diff_kmers(reads, varargin)

global organism;







chromosomes = organism.Chromosomes;
transcripts = organism.Transcripts;
exons = organism.Exons;

normalization = 'none';
read_len = 50;
granularity = 2000;

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'ReadLen')
		read_len = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

if isnan(read_len)
	error 'Read length must be specified.';
end

S = length(reads.Raw);

expr = struct;
expr.Chromosomes = cell(length(chromosomes.Name), S);

expr.Meta = reads.Meta;
expr.Meta.Type = 'Exon expression';
expr.Meta.Organism = organism.Name;
expr.Meta.TotalSeqReads = zeros(S, 1);

for s = 1:S
	sample_id = reads.Meta.Sample.ID{s};
	
	fprintf(1, ['Constructing genomic expression track based on RNA-seq ' ...
		'sample %s...\n'], sample_id);
	al = align_reads(filter_query(reads, s), 'genome', ...
		'MaxMismatches', 2, 'AllowAlignments', 1, ...
		'Columns', 'target,offset', varargin{:});
	
	chroms = chromosome_sym2num(al.Target);
	offsets = al.Offset + round(read_len / 2);
	
	for chr = 1:length(chromosomes.Name)
		expr.Chromosomes{chr, s} = histc(offsets(chroms == chr), ...
			0:granularity:chromosomes.Length(chr)+1);
		
		% Enable this code if you want overlapping windows.
		%odd_track = histc(offsets(chroms == chr), ...
		%	0:granularity:chromosomes.Length(chr)+1);
		%even_track = histc(offsets(chroms == chr), ...
		%	granularity/2:granularity:chromosomes.Length(chr)+1-granularity/2);
		
		%combined_track = zeros(1, length(even_track) + length(odd_track));
		%combined_track(2*(0:length(even_track)-1)+1) = even_track;
		%combined_track(2*(1:length(even_track)) = odd_track;
		%expr.Chromosomes{chr, s} = 
		
	end
	
	expr.Meta.TotalSeqReads(s) = sum(al.TotalReads);
end



