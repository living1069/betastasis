
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
%    MIRNA_EXPRESSION_RNASEQ(..., 'Trim', LEN) tells the function to trim all
%    reads to a length of LEN nucleotides. Nucleotides are trimmed from the %    end of the reads. Defaults to trimming reads to 18 bases.
%
%    MIRNA_EXPRESSION_RNASEQ(..., 'LengthFilter', [MIN_LEN]) tells the function
%    to filter out all reads shorter than MIN_LEN before alignment. Defaults
%    to no filtering.
%
%    MIRNA_EXPRESSION_RNASEQ(..., 'LengthFilter', [MIN MAX]) tells the function
%    to filter out all reads outside the length range [MIN MAX]. Defaults to
%    no filtering.
%
%    MIRNA_EXPRESSION_RNASEQ(..., 'Normalization', NORM) tells the function
%    to normalize the produced miRNA expression values using the normalization
%    method NORM. Valid selections include 'rpkm' and 'none' (the default).
% 

function expr = mirna_expression_rnaseq(reads, varargin)

global organism;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

max_mismatches = 0;
normalization = 'none';
trim_len = 18;
filter_len = [];

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MaxMismatches')
		max_mismatches = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Trim')
		trim_len = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'LengthFilter')
		filter_len = varargin{k+1};
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

seq_files = seq_resource_files(reads);

% Figure out what kind of data we're dealing with.
[color, quality] = seq_read_type(seq_files{1});
if quality, error 'FASTQ files not supported yet.'; end

pre_mirna_name_to_idx = containers.Map(pre_mirnas.Name, ...
	num2cell(1:length(pre_mirnas.Name)));

expr = struct;
expr.miRNA = zeros(length(mirnas.Name), length(seq_files));
expr.pre_miRNA = zeros(length(pre_mirnas.Name), length(seq_files));

expr.Meta = struct;
if isstruct(reads), expr.Meta = reads.Meta; end

expr.Meta.Type = 'miRNA expression';
expr.Meta.Organism = organism.Name;
expr.Meta.miRNAVersion = organism.miRNAVersion;
expr.Meta.Normalization = repmat({normalization}, length(seq_files), 1);
expr.Meta.MismatchesAllowed = ones(length(seq_files), 1) * max_mismatches;
expr.Meta.TotalSeqReads = zeros(length(seq_files), 1);

for seq_file = 1:length(seq_files)
	fprintf(1, 'Aligning reads to miRBase using Bowtie...\n');
	al = align_reads(seq_files{seq_file}, 'pre_mirnas', ...
		'MaxMismatches', max_mismatches, 'AllowAlignments', 10, ...
		'Trim', trim_len, 'LengthFilter', filter_len, ...
		'Columns', 'read,ref,offset,sequence');
	
	expr.Meta.TotalSeqReads(seq_file) = al.TotalReads;

	read_ids = al.ReadID;
	match_names = al.Target;
	read_offsets = al.Offset;
	read_sequences = al.Sequence;

	read_lengths = zeros(length(read_sequences), 1);
	for k = 1:length(read_sequences)
		read_lengths(k) = length(read_sequences{k});
	end
	clear read_sequences;

	mirna_lengths = zeros(length(mirnas.Sequence), 1);
	for k = 1:length(mirna_lengths)
		mirna_lengths(k) = length(mirnas.Sequence{k});
	end

	fprintf(1, 'Summarizing miRNA expression levels...\n');

	for k = 1:length(match_names)
		p = pre_mirna_name_to_idx(match_names{k});
		
		% Check if the read is aligned to a spot within one of the mature miRNA
		% sequences. If it is, we assume that the read comes from the mature
		% miRNA and not from a pre-miRNA molecule. In reality, the read
		% could very well originate from the pre-miRNA, but there's no way
		% to tell.
		found_mirna_match = 0;
		for m = 1:pre_mirnas.MatureCount(p)
			idx = pre_mirnas.Matures(p, m);
			if read_offsets(k) + 1 >= pre_mirnas.MatureOffsets(p, m) && ...
			   read_offsets(k) + 1 + read_lengths(k) <= ...
			   pre_mirnas.MatureOffsets(p, m) + mirna_lengths(idx)
				expr.miRNA(idx, seq_file) = expr.miRNA(idx, seq_file) + 1;
				found_mirna_match = 1;
				break;
			end
		end
			
		if found_mirna_match == 0
			expr.pre_miRNA(p, seq_file) = expr.pre_miRNA(p, seq_file) + 1;
		end
	end
	
	if strcmp('none', normalization)
		continue;
	end

	fprintf(1, 'Performing RPKM normalization on pre-miRNA expression levels...\n');
	pre_mirna_kbp = zeros(length(pre_mirnas.Sequence), 1);
	for k = 1:length(pre_mirnas.Sequence)
		pre_mirna_kbp(k) = length(pre_mirnas.Sequence{k}) / 1000;
	end

	expr.pre_miRNA(:, seq_file) = expr.pre_miRNA(:, seq_file) ./ ...
		pre_mirna_kbp / (reads_processed / 1e6);

	% We don't normalize mature miRNA by transcript length, since their
	% lengths are so uniform.
	fprintf(1, ['Performing RPKM normalization on mature miRNA expression ' ...
				'levels...\n']);
	expr.miRNA(:, seq_file) = expr.miRNA(:, seq_file) / (reads_processed / 1e6);
end

