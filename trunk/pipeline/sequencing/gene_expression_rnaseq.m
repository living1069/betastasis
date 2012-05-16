
% GENE_EXPRESSION_RNASEQ    Calculate gene expression from RNA-seq reads.
%
%    EXPR = GENE_EXPRESSION_RNASEQ(READS) takes as input a realized dataset
%    READS containing short sequence reads. The function outputs a gene
%    expression profile for each sample in READS.
%
%    GENE_EXPRESSION_RNASEQ(..., 'Normalization', NORM) tells the function
%    to normalize the produced gene expression values using the normalization
%    method NORM. Valid selections include 'rpkm' and 'none' (the default).
%
%    All expression values are returned in the natural scale and are not log-
%    transformed. The rows of the expression matrix in EXPR correspond
%    one-to-one with gene names in organism.Genes.Name.
%
%    Also see HELP ALIGN_READS for a number of additional options for the
%    alignment step.

function expr = gene_expression_rnaseq(reads, varargin)

global organism;
genome = organism.Genes;
transcriptome = organism.Transcripts;

S = length(reads.Raw);

expr = struct;
expr.Mean = zeros(length(genome.Name), S);

expr.Meta = struct;
if isstruct(reads), expr.Meta = reads.Meta; end

expr.Meta.Type = 'Gene expression';
expr.Meta.Organism = organism.Name;
%expr.Meta.OrganismVersion = organism.Version;
expr.Meta.TotalSeqReads = zeros(S, 1);


for s = 1:S
	if isfield(reads.Meta.Sample, 'ID')
		sample_id = reads.Meta.Sample.ID{s};
	else
		sample_id = reads.Meta.Sample.Filename{s};
	end
	
	fprintf(1, 'Calculating gene expressions for RNA-seq sample %s...\n', ...
		sample_id);
	al = align_reads(filter_query(reads, s), 'transcripts', ...
		'MaxMismatches', 2, 'AllowAlignments', 20, ...
		'Columns', 'read,target', varargin{:});
	
	read_ids = al.ReadID;
	transcripts = al.Target;
	
	transcript_expr = zeros(length(transcriptome.Name), 1);
	transcript_indices = transcript_idx(transcripts);

	fprintf(1, 'Summarizing transcript expression levels...\n');
	run_ends = [ find(strcmp(read_ids(1:end-1), read_ids(2:end)) == 0); ...
		length(read_ids) ];
	run_lengths = diff([0; run_ends]);

	pos = 1;
	for r = 1:length(run_lengths)
		for k = pos:pos+run_lengths(r)-1
			t = transcript_indices(k);
			transcript_expr(t) = transcript_expr(t) + 1 / run_lengths(r);
		end
		pos = pos + run_lengths(r);
	end
	
	fprintf('Summarizing gene expression from transcript expression...\n');
	transcript_expr = struct('Mean', transcript_expr);
	gene_expr = gene_expression_from_transcript_expression( ...
		transcript_expr, 'sum');
	expr.Mean(:, s) = gene_expr.Mean;
	expr.Meta.TotalSeqReads(s) = sum(al.TotalReads);
end

