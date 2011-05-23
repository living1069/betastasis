
% MIRNA_EXPRESSION_RNASEQ    Calculate microRNA expression from RNA-seq reads.
%
%    EXPR = MIRNA_EXPRESSION_RNASEQ(READS) takes as input a realized dataset
%    READS containing sequencer reads. The function outputs a microRNA
%    expression dataset that contains both pre-miRNA and mature miRNA.
%
%    MIRNA_EXPRESSION_RNASEQ(..., 'MaxMismatches', MM) sets the maximum number
%    of mismatches allowed while aligning reads against microRNA transcripts.
%    Defaults to 0 mismatches.
%
%    MIRNA_EXPRESSION_RNASEQ(..., 'Normalization', NORM) tells the function
%    to normalize the produced miRNA expression values using the normalization
%    method NORM. Valid selections include 'rpkm' and 'none' (the default).

% Author: Matti Annala <matti.annala@tut.fi>

function [mature_expr, pre_expr] = mirna_expression_rnaseq(reads, varargin)

global organism;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

max_mismatches = 0;
normalization = 'none';
trim_len = 18;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MaxMismatches')
		max_mismatches = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Normalization')
		if strcmpi('none', varargin{k+1})
			normalization = 'none';
		elseif strcmpi('rpkm', varargin{k+1})
			normalization = 'rpkm';
		else
			error(['Invalid normalization method specified. Valid ' ...
				   'alternatives include "none" and "rpkm".']);
		end
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

S = length(reads.Raw);

pre_mirna_name_to_idx = containers.Map(pre_mirnas.Name, ...
	num2cell(1:length(pre_mirnas.Name)));
	
mirna_lengths = zeros(length(mirnas.Sequence), 1);
for k = 1:length(mirna_lengths)
	mirna_lengths(k) = length(mirnas.Sequence{k});
end

mature_expr.Mean = zeros(length(mirnas.Name), S);
pre_expr.Mean = zeros(length(pre_mirnas.Name), S);

mature_expr.Meta = reads.Meta;
mature_expr.Meta.Type = 'miRNA expression';
mature_expr.Meta.Organism = organism.Name;
mature_expr.Meta.DateCreated = repmat({datestr(now)}, S, 1);
mature_expr.Meta.Normalization = repmat({normalization}, S, 1);
mature_expr.Meta.MismatchesAllowed = ones(S, 1) * max_mismatches;
mature_expr.Meta.TotalSeqReads = zeros(S, 1);

pre_expr.Meta = mature_expr.Meta;
pre_expr.Meta.Type = 'pre-miRNA expression';

for s = 1:length(reads.Raw)
	fprintf(1, 'Calculating microRNA expression levels for sample %s:\n', ...
		reads.Meta.Sample.Filename{s});
	
	extracted = extract_reads(filter_query(reads, s));
	
	trimmed = trim_reads(extracted, trim_len);
	
	fprintf(1, '-> Aligning reads to miRBase using Bowtie...\n');
	al = align_reads(trimmed, 'pre_mirnas', ...
		'MaxMismatches', max_mismatches, 'AllowAlignments', 10, ...
		'Columns', 'ref,offset,sequence');
	
	mature_expr.Meta.TotalSeqReads(s) = al.TotalReads;

	match_names = al.Target;
	read_offsets = al.Offset;

	read_lengths = zeros(length(al.Sequence), 1);
	for k = 1:length(al.Sequence)
		read_lengths(k) = length(al.Sequence{k});
	end

	fprintf(1, '-> Summarizing miRNA expression levels...\n');
	
	pre_idx = cell2mat(pre_mirna_name_to_idx.values(match_names));
	
	for k = 1:length(match_names)
		p = pre_idx(k);
		
		% Check if the read is aligned to a spot within one of the mature miRNA
		% sequences. If it is, we assume that the read comes from the mature
		% miRNA and not from a pre-miRNA molecule. In reality, the read
		% could very well originate from the pre-miRNA, but there's no way
		% to tell.
		found_mirna_match = 0;
		for m = 1:pre_mirnas.MatureCount(p)
			idx = pre_mirnas.Matures(p, m);
			if read_offsets(k) >= pre_mirnas.MatureOffsets(p, m) && ...
			   read_offsets(k) + read_lengths(k) <= ...
			   pre_mirnas.MatureOffsets(p, m) + mirna_lengths(idx)
				mature_expr.Mean(idx, s) = mature_expr.Mean(idx, s) + 1;
				found_mirna_match = 1;
				break;
			end
		end
			
		if found_mirna_match == 0
			pre_expr.Mean(p, s) = pre_expr.Mean(p, s) + 1;
		end
	end
	
	if strcmp('none', normalization), continue, end

	fprintf(1, '-> Normalizing pre-miRNA expression levels using RPKM...\n');
	pre_expr.Mean(:, s) = pre_expr.Mean(:, s) / (reads_processed / 1e6);

	% We don't normalize mature miRNA by transcript length, since their
	% lengths are so uniform.
	fprintf(1, '-> Normalizing mature miRNA expression levels using RPKM...\n');
	mature_expr.Mean(:, s) = mature_expr.Mean(:, s) / (reads_processed / 1e6);
end

