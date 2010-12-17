
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

normalization = 'none';

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Normalization')
		if strcmpi('none', varargin{k+1})
			normalization = 'None';
		elseif strcmpi('rpkm', varargin{k+1})
			normalization = 'RPKM';
		elseif strcmpi('quantile', varargin{k+1})
			normalization = 'Quantile';
		elseif regexpi(varargin{k+1}, 'rpkm.*quantile')
			normalization = 'RPKM + Quantile';
		else
			error 'Unrecognized normalization option given as a parameter.';
		end
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

seq_files = seq_resource_files(reads);

expr = struct;
expr.Mean = zeros(length(genome.Name), length(seq_files));

expr.Meta = struct;
if isstruct(reads), expr.Meta = reads.Meta; end

expr.Meta.Type = 'Gene expression';
expr.Meta.Organism = organism.Name;
expr.Meta.OrganismVersion = organism.Version;
expr.Meta.Normalization = repmat({ normalization }, length(seq_files), 1);
expr.Meta.TotalSeqReads = zeros(length(seq_files), 1);

for seq_file = 1:length(seq_files)
	if isfield(reads.Meta.Sample, 'ID')
		sample_id = reads.Meta.Sample.ID{seq_file};
	else
		sample_id = reads.Meta.Sample.Filename{seq_file};
	end
	
	fprintf(1, 'Calculating gene expressions for RNA-seq sample %s...\n', ...
		sample_id);
	al = align_reads(seq_files{seq_file}, 'transcripts', ...
		'MaxMismatches', 2, 'AllowAlignments', 20, ...
		'Columns', 'read,target', varargin{:});
	
	read_ids = al.ReadID;
	transcripts = al.Target;
	
	transcript_expr = zeros(length(transcriptome.Name), 1);

	fprintf(1, 'Summarizing transcript expression levels...\n');
	run_ends = [ find(strcmp(read_ids(1:end-1), read_ids(2:end)) == 0); ...
		length(read_ids) ];
	run_lengths = diff([0; run_ends]);

	transcript_indices = transcript_idx(transcripts);

	pos = 1;
	for r = 1:length(run_lengths)
		for k = pos:pos+run_lengths(r)-1
			t = transcript_indices(k);
			transcript_expr(t) = transcript_expr(t) + 1 / run_lengths(r);
		end
		pos = pos + run_lengths(r);
	end

	if regexp(normalization, 'RPKM')
		fprintf(1, ['Performing RPKM normalization on transcript ' ...
		            'expression levels...\n']);
		
		transcript_kbp = zeros(length(transcriptome.Sequence), 1);
		for k = 1:length(transcriptome.Sequence)
			transcript_kbp(k) = length(transcriptome.Sequence{k}) / 1000;
		end

		transcript_expr = transcript_expr ./ transcript_kbp / ...
			(al.TotalReads / 1e6);
	end
	
	fprintf(1, ['Summarizing a gene expression profile from transcript ' ...
	            'expression levels...\n']);
	transcript_expr = struct('Mean', transcript_expr);
	gene_expr = gene_expression_from_transcript_expression( ...
		transcript_expr, 'sum');
	expr.Mean(:, seq_file) = gene_expr.Mean;
	expr.Meta.TotalSeqReads(seq_file) = al.TotalReads;
end

if regexp(normalization, 'Quantile')
	fprintf(1, ['Performing quantile normalization on gene ' ...
	            'expression levels...\n']);
	expr.Mean = quantilenorm(expr.Mean);
end

