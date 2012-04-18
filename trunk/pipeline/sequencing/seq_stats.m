function [] = seq_stats(reads, varargin)

global organism;
transcripts = organism.Transcripts;

report_dir = '~/rnaseq_stats';
[~, ~] = mkdir(report_dir);

tx_len = zeros(length(transcripts.Sequence), 1);
for k = 1:length(transcripts.Sequence)
	tx_len(k) = length(transcripts.Sequence{k});
end

rrna_genes = { 'LOC1000008587'; 'LOC1000008588'; 'LOC1000008589' };
rrna_genes = cat(1, rrna_genes, ...
	organism.Genes.Name(rx(organism.Genes.Name, 'RN5S')));

rrna_transcripts = [];
for g = rrna_genes'
	rrna_transcripts(end+1) = organism.Genes.Transcripts(g, ...
		1:organism.Genes.TranscriptCount);
end

extracted = extract_reads(reads);

%for s = 1:length(extracted.Raw)
for s = 1
	
	% First we align against the transcriptome.
	[al, tx_unaligned] = align_reads(filter_query(extracted, s), ...
		'transcripts', 'Columns', 'target,offset,sequence', ...
		'ReportAlignments', 20, 'AllowAlignments', 20, varargin{:});
	
	fprintf(1, 'Total reads aligned to transcriptome: %d\n', al.
	
	target = transcript_idx(al.Target);
	read_len = zeros(length(target), 1);
	for k = 1:length(al.Sequence), read_len(k) = length(al.Sequence{k}); end
	
	normalized_offsets = offset ./ (tx_len(target) - read_len + 1);
	figure; hist(normalized_offsets, 0:0.01:1); xlim([0 1]);
	saveas(gcf, sprintf('%s/%d_norm_offsets.pdf', report_dir, seq_file));
	fprintf(1, 'Normalized read offset distribution: %s\n', ...
		sprintf('%s/%d_norm_offsets.pdf', report_dir, seq_file));
	
	figure; hist(offset, 1:max(offset)); xlim([0 max(offset)]);
	saveas(gcf, sprintf('%s/%d_5p_offsets.pdf', report_dir, seq_file));
	fprintf(1, '5'' read offset distribution: %s\n', ...
		sprintf('%s/%d_5p_offsets.pdf', report_dir, seq_file));
	
	offsets_3p = tx_len(target) - (offset + read_len) + 2;
	figure; hist(offsets_3p, 1:max(offsets_3p)); xlim([0 max(offsets_3p)]);
	saveas(gcf, sprintf('%s/%d_3p_offsets.pdf', report_dir, seq_file));
	fprintf(1, '3'' read offset distribution: %s\n', ...
		sprintf('%s/%d_3p_offsets.pdf', report_dir, seq_file));
	
	cds_left = organism.Transcripts.CDS(target, 1);
	cds_right = organism.Transcripts.CDS(target, 2);
	in_cds = (offset >= cds_left) & (offset + read_len - 1 <= cds_right);
	fprintf(1, 'Percentage of reads within CDS: %.1f\n', ...
		sum(in_cds) / length(in_cds) * 100);
	
	total_exon_len = 0;
	total_cds_len = 0;
	
	for m = 1:length(target)
		idx = tx_list(m);
		total_exon_len = total_exon_len + ...
			length(organism.Transcripts.Sequence{idx});
		if ~any(isnan(organism.Transcripts.CDS(idx, :)))
			total_cds_len = total_cds_len + ...
				organism.Transcripts.CDS(idx, 2) - ...
				organism.Transcripts.CDS(idx, 1) + 1;
		end
	end
			
	fprintf(1, 'Expected percentage of reads within CDS: %.1f\n', ...
		total_cds_len / total_exon_len * 100);

		
		
		
		
	% Then we align the remaining reads against the genome.
	al = align_reads(tx_unaligned, 'genome', ...
		'ReportAlignments', 10, 'AllowAlignments', 10, ...
		'Columns', 'target,offset,sequence', varargin{:});

end

