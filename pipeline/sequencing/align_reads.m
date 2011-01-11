
% ALIGN_READS       Align sequence reads against reference sequences
%
%    ALIGNMENTS = ALIGN_READS(READS, INDEX, ...) aligns the reads given in the
%    input file READS against the reference sequences specified by INDEX. The
%    input file format must be one of FASTA, FASTQ, CSFASTA, CSFASTQ or SMS.
%    Valid values for INDEX include 'genome', 'transcripts', 'exons', 'mirnas'
%    and 'pre_mirnas', which all refer to the corresponding annotations for the
%    currently selected organism (see HELP ORGANISM). Alternatively, the index
%    may be specified as a cell array of sequences, or a structure with 
%    fields 'Name' and 'Sequence'.
%    
%    ALIGN_READS(..., 'Aligner', ALIGNER) specifies the alignment algorithm
%    that will be used for aligning the reads against the reference sequences.
%    The default is to use Helisphere for SMS files, and Bowtie for everything
%    else.
%
%    ALIGN_READS(..., 'MaxMismatches', MM) tells the alignment algorithm
%    to only include matches with MM mismatches or less. The number of
%    mismatches is specified in nucleotide space or colorspace, depending on
%    the data type used.
%
%    ALIGN_READS(..., 'Trim', TRIM_LEN) tells the alignment algorithm to
%    first trim all reads so that only the first TRIM_LEN nucleotides
%    (or colors) from the 5' end are aligned.
%
%    ALIGN_READS(..., 'LengthFilter', LEN_RANGE) tells the alignment algorithm
%    to discard all reads outside the range specified by the two element vector
%    LEN_RANGE.
%
%    ALIGN_READS(..., 'MaxThreads', MAX_THR) constrains the alignment algorithm
%    to use a maximum of MAX_THR processing threads. The default is to use 
%    pipeline_config.MaxThreads threads. Not all alignment algorithms
%    support multithreading; for those algorithms this option has no effect.

function alignments = align_reads(reads, index, varargin)

aligner = '';

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Aligner')
		aligner = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

if isempty(aligner)
	if regexpi(reads, '.*\.sms')
		alignments = helisphere_align(reads, index, varargin{:});
	else
		alignments = bowtie_align2(reads, index, varargin{:});
	end
elseif strcmpi(aligner, 'bowtie') && regexpi(reads, '.*\.sms')
	fprintf(1, 'Converting SMS reads to FASTA...\n');
	fasta_tmp = ptemp;
	[status, ~] = unix(sprintf(['%s/tools/helisphere/bin/sms2txt ' ...
		'--fasta --input_file %s --output_file_prefix %s'], ...
		ppath, reads, fasta_tmp));
	if status ~= 0, error 'SMS to FASTA conversion failed.'; end
	
	fasta_tmp = [fasta_tmp '.fa'];  % sms2txt appends .fa to the filename
	alignments = bowtie_align2(fasta_tmp, index, varargin{:});
	delete(fasta_tmp);
elseif strcmpi(aligner, 'bowtie')
	alignments = bowtie_align2(reads, index, varargin{:});
elseif strcmpi(aligner, 'helisphere')
	alignments = helisphere_align(reads, index, varargin{:});
elseif strcmpi(aligner, 'gassst')
	alignments = gassst_align(reads, index, varargin{:});
else
	error 'Requested alignment software is not supported.';
end

