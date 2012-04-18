
% CHIPSEQ_PEAKS     Call protein binding sites based on ChIP-seq reads
%
%    PEAKS = CHIPSEQ_PEAKS(TEST, REF) makes protein binding site calls based on
%    the ChIP-seq experiments TEST and REF, where the protein of interest was
%    crosslinked in TEST, and REF is a control experiment with a non-specific
%    protein (or no protein).
%
%    CHIPSEQ_PEAKS(..., 'PromoterRegion', REGION) can be used to adjust the
%    expected position of a gene's promoter region. REGION must be specified
%    as a two element vector that specifies the position of a gene's promoter
%    region relative to its transcription start site (TSS). The range is
%    specified inclusively at both ends. The default value for this option
%    is [-5000, 0].
%
%    CHIPSEQ_PEAKS(..., 'PeakCaller', CALLER) uses the peak calling algorithm
%    CALLER to determine binding sites. The function defaults to using 'MACS'
%    (which is currently the only supported algorithm).

% Author: Matti Annala <matti.annala@tut.fi>

function peaks = chipseq_peaks(test, ref, varargin)

global organism;
genes = organism.Genes;

peak_caller = 'MACS';
promoter_region = [-5000 0];

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'PeakCaller')
		peak_caller = varargin{k+1};
		if strcmpi(peak_caller, 'MACS')
			error 'Only MACS is currently supported as a peak caller.';
		end
		drop_args(k:k+1) = true;
		continue;
	end
	
	if strcmpi(varargin{k}, 'PromoterRegion')
		promoter_region = varargin{k+1};
		if numel(promoter_region) ~= 2
			error 'Promoter region must be specified as a relative interval.';
		end
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

test_files = seq_resource_files(test)
ref_files = seq_resource_files(ref)

if length(test_files) ~= length(ref_files)
	error 'Equally many test and reference samples must be provided.';
end

peaks = struct;
peaks.Meta = struct;

if isstruct(reads), peaks.Meta = reads.Meta; end

expr.Meta.Type = 'ChIP-seq peaks';
expr.Meta.TotalSeqReads = zeros(length(seq_files), 1);

for seq_file = 1:length(test_files)
	tmp = ptemp;
	
	% Calculate the average read length (tag size).
	test_fasta = FastaFile(test_files{seq_file});
	ref_fasta = FastaFile(ref_files{seq_File});
	
	if abs(test_fasta.AvgReadLength - ref_fasta.AvgReadLength) > 5
		fprintf(1, ['WARNING: Large difference in read lengths for the ' ...
			'test and reference samples. Trimming is adviced...\n']);
	end
	
	avg_tag_size = round(mean( ...
		[test_fasta.AvgReadLength, ref_fasta.AvgReadLength]));
	
	fprintf(1, 'Aligning test sample reads against the genome...\n');
	test_al = align_reads(test_files{seq_file}, 'genome', ...
		'MaxMismatches', 2, 'AllowAlignments', 10, 'ReportAlignments', 10, ...
		'Columns', 'offset,read,target', varargin{:});
		
	fprintf(1, 'Aligning reference sample reads against the genome...\n');
	ref_al = align_reads(ref_files{seq_file}, 'genome', ...
		'MaxMismatches', 2, 'AllowAlignments', 10, 'ReportAlignments', 10, ...
		'Columns', 'offset,read,target', varargin{:});

	fprintf(1, 'Writing alignments to disk in BED format...\n');
	test_bed = fopen([tmp '.test.bed'], 'W');
	for k = 1:length(test_al.Offset)
		fprintf(test_bed, '%s\t%d\t%d\n', ...
			organism.Chromosomes.Name{test_al.Chromosome}, ...
			test_al.Offset(k), test_al.Offset(k) + length(test_al.Sequence{k}));
	end
	fclose(test_bed);
	
	ref_bed = fopen([tmp '.ref.bed'], 'W');
	for k = 1:length(ref_al.Offset)
		fprintf(ref_bed, '%s\t%d\t%d\n', ...
			organism.Chromosomes.Name{ref_al.Chromosome}, ...
			ref_al.Offset(k), ref_al.Offset(k) + length(ref_al.Sequence{k}));
	end
	fclose(ref_bed);

	fprintf(1, 'Finding ChIP-seq peaks based on alignments...\n');
	[status, ~] = unix(sprintf( ...
		['PYTHONPATH=%s/tools/macs/lib %s/tools/macs/bin/macs ' ...
		 '--format=BED -tsize=%d -t %s -c %s'], ...
		ppath, ppath, avg_tag_size, test_bed, ref_bed));
		
	safe_delete([tmp '.test.bed']);
	safe_delete([tmp '.ref.bed']);
	
	transcript_indices = transcript_idx(transcripts);

end

