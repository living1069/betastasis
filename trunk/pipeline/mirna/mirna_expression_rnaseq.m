
% MIRNA_EXPRESSION_RNASEQ    Calculate microRNA expression from RNA-seq reads.
%
%    EXPR = MIRNA_EXPRESSION_RNASEQ(READS) takes as input a realized dataset
%    READS containing sequencer reads. The function outputs a microRNA
%    expression dataset that contains both pre-miRNA and mature miRNA.
%
%    MIRNA_EXPRESSION_RNASEQ(..., 'MaxMismatches', MM) sets the maximum number
%    of mismatches allowed while aligning reads against microRNA transcripts.
%    Defaults to 0 mismatches.

% Author: Matti Annala <matti.annala@tut.fi>

function mirna_expr = mirna_expression_rnaseq(reads, varargin)

global organism;
mirnas = organism.miRNA;

max_mismatches = 0;
trim_len = 18;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'max.*mismatch')
		max_mismatches = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

S = length(reads.url);

mirna_name_to_idx = containers.Map(mirnas.Name, ...
	num2cell(1:length(mirnas.Name)));
	
mirna_expr.mean = zeros(length(mirnas.Name), S);

mirna_expr.meta = reads.meta;
mirna_expr.meta.type = 'miRNA expression';
mirna_expr.meta.organism = organism.Name;
mirna_expr.meta.mismatches_allowed = ones(1, S) * max_mismatches;
mirna_expr.meta.total_reads = nan(1, S);

for s = 1:S
	fprintf('Calculating microRNA expression levels for sample %s:\n', ...
		reads.meta.sample_id{s});
	
	trimmed = trim_reads(filter(reads, s), trim_len);
	
	fprintf('-> Aligning reads to miRBase using Bowtie...\n');
	alignments = bowtie_align(trimmed, 'mirnas', ...
		sprintf('-v%d -m10', max_mismatches));
	
	mirna_expr.meta.total_reads(s) = alignments.total_reads;
	
	fprintf('-> Summarizing miRNA expression levels...\n');
	for al = iterate_alignments(alignments)
		mirna_idx = cell2mat(mirna_name_to_idx.values(al.target));
		for p = mirna_idx'
			mirna_expr.mean(p, s) = mirna_expr.mean(p, s) + 1;
		end
	end
end

