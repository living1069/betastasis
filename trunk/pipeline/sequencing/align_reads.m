
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
%    ALIGN_READS(..., 'MaxThreads', MAX_THR) constrains the alignment algorithm
%    to use a maximum of MAX_THR processing threads. The default is to use 
%    pipeline_config.MaxThreads threads. Not all alignment algorithms
%    support multithreading; for those algorithms this option has no effect.

function [alignments, unaligned] = align_reads(reads, index, varargin)

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

if ~isfield(reads, 'Raw')
	error 'Reads must be passed as a dataset.';
end
if length(reads.Raw) ~= 1
	error 'align_reads() can only be called one sample at a time.';
end

if isempty(aligner)
	if regexpi(reads.Meta.Sequence.Format{1}, 'SMS')
		align = @helisphere_align;
	else
		align = @bowtie_align2;
	end
elseif strcmpi(aligner, 'bowtie')
	align = @bowtie_align2;
elseif strcmpi(aligner, 'helisphere')
	align = @helisphere_align;
elseif strcmpi(aligner, 'gassst')
	align = @gassst_align;
else
	error 'Requested alignment software is not supported.';
end

if nargout == 2
	[alignments, unaligned] = align(reads, index, varargin{:});
elseif nargout == 1
	alignments = align(reads, index, varargin{:});
end

