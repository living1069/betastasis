
% NORMALIZE_EXPR    Normalize gene expression values across samples
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

function expr = normalize_expr(expr, method)

global organism;
genome = organism.Genes;
transcriptome = organism.Transcripts;

if nargin < 2, method = 'none'; end

S = size(expr.Mean, 2);

expr = struct;
expr.Mean = zeros(length(genome.Name), S);

expr.Meta = struct;
if isstruct(reads), expr.Meta = reads.Meta; end

expr.Meta.Type = 'Gene expression';
expr.Meta.Organism = organism.Name;
%expr.Meta.OrganismVersion = organism.Version;
expr.Meta.Normalization = repmat({ normalization }, S, 1);
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
	
	if regexpi(normalization, 'RPKM')
		if ~isnan(al.TotalReads)
			fprintf(1, ['Performing RPKM normalization on transcript ' ...
						'expression levels...\n']);
			
			transcript_kbp = zeros(length(transcriptome.Sequence), 1);
			for k = 1:length(transcriptome.Sequence)
				transcript_kbp(k) = length(transcriptome.Sequence{k}) / 1000;
			end
			
			transcript_expr = transcript_expr ./ transcript_kbp / ...
				(al.TotalReads / 1e6);
		else
			fprintf(1, ['Cannot perform RPKM normalization. No total read ' ...
				'count information available.\n']);
		end
	end
	
	fprintf(1, ['Summarizing a gene expression profile from transcript ' ...
	            'expression levels...\n']);
	transcript_expr = struct('Mean', transcript_expr);
	gene_expr = gene_expression_from_transcript_expression( ...
		transcript_expr, 'sum');
	expr.Mean(:, s) = gene_expr.Mean;
	expr.Meta.TotalSeqReads(s) = sum(al.TotalReads);
end

