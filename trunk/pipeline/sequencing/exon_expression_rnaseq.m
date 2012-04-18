
function expr = exon_expression_rnaseq(reads, varargin)

global organism;
transcripts = organism.Transcripts;
exons = organism.Exons;

normalization = 'none';
read_len = nan;

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

S = length(reads.url);

expr = struct;
expr.mean = zeros(length(exons.ID), S);
expr.total_seq_reads = nan(S, 1);
expr.meta = reads.meta;
expr.meta.type = 'Exon expression';
expr.meta.organism = organism.Name;


for s = 1:S
	fprintf(1, 'Calculating exon expressions for RNA-seq sample %s...\n', ...
		reads.meta.sample_id{s});
	alignments = bowtie_align(filter(reads, s), 'transcripts', '-v2 -m20');
	
	al = all_alignments(alignments);
	
	targets = transcript_idx(al.target);
	offsets = al.offset;
	
	exon_expr = zeros(length(exons.ID), 1);

	fprintf(1, 'Summarizing exon expression levels...\n');
	for k = 1:length(targets)
		exon_pos = transcripts.ExonPos{targets(k)};
		if isempty(exon_pos), continue, end
			
		read_pos = [offsets(k), offsets(k) + read_len - 1];
		overlap_exons = transcripts.Exons{targets(k)}( ...
			read_pos(1) <= exon_pos(:, 2) - 5 & ...
			read_pos(2) >= exon_pos(:, 1) + 5);
		
		exon_expr(overlap_exons) = exon_expr(overlap_exons) + 1;
	end
	
	expr.mean(:, s) = exon_expr;
	expr.total_seq_reads(s) = sum(alignments.total_reads);
end

