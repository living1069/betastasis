function [] = rnaseq_stats(reads, report_dir, varargin)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;

read_len = 90;

contaminants = human_contaminants();
contaminant_name_to_idx = containers.Map(contaminants.Name, ...
	num2cell(1:length(contaminants.Name)));
	
[~, ~] = mkdir(report_dir);
	
S = length(reads.Raw);

if isfield(reads.Meta.Sample, 'ID')
	sample_ids = reads.Meta.Sample.ID;
else
	sample_ids = reads.Meta.Sample.Filename;
end

total_reads = zeros(1, S);
contaminant_reads = zeros(length(contaminants.Name), S);
exonic_reads = zeros(1, S);
intronic_reads = zeros(1, S);
intergenic_reads = zeros(1, S);
unknown_reads = zeros(1, S);

for s = 1:S
	fprintf(1, 'Checking RNA contamination in sample %s...\n', sample_ids{s});
	[al, unaligned] = align_reads(filter_query(reads, s), contaminants, ...
		'MaxMismatches', 2, 'ReportAlignments', 1, ...
		'Columns', 'target', varargin{:});
		
	total_reads(s) = sum(al.TotalReads);
	
	target = cell2mat(contaminant_name_to_idx.values(al.Target));
	for t = target'
		contaminant_reads(t, s) = contaminant_reads(t, s) + 1;
	end
	
	fprintf(1, 'Checking exonic reads in sample %s...\n', sample_ids{s});
	[al, unaligned] = align_reads(unaligned, 'transcripts', ...
		'MaxMismatches', 2, 'ReportAlignments', 1, ...
		'Columns', 'target,offset', varargin{:});

	tx_len = nan(length(transcripts.Sequence), 1);
	for k = 1:length(transcripts.Sequence)
		tx_len(k) = length(transcripts.Sequence{k});
	end

	target = transcript_idx(al.Target);
	offset = double(al.Offset);
	
	exonic_reads(s) = length(target);
	
	normalized_offsets = offset ./ (tx_len(target) - read_len + 1);
	figure; hist(normalized_offsets, 0:0.01:1); xlim([0 1]);
	saveas(gcf, sprintf('%s/%d_norm_offsets.pdf', report_dir, s));
	
	tx_len_alignments = hist(tx_len(target), 1:100:max(tx_len));
	figure; bar(log10(tx_len_alignments / sum(tx_len_alignments)));
	xlabel('Transcript length'); ylabel('Percent aligned reads');
	saveas(gcf, sprintf('%s/%d_tx_len_alignments.pdf', report_dir, s));
	
	fprintf(1, 'Checking genomic reads in sample %s...\n', sample_ids{s});
	al = align_reads(unaligned, 'genome', ...
		'MaxMismatches', 2, 'ReportAlignments', 1, ...
		'Columns', 'target,offset', varargin{:});

	chr = chromosome_sym2num(al.Target);
	offset = al.Offset;
	
	gene_pos = organism.Genes.Position;
	gene_chr = organism.Genes.Chromosome;
	
	% Check if the read is intragenic. If it is, we count it as an intronic
	% read, since exonic reads were already discovered during the previous step.
	for c = 1:max(chr)
		chr_genes = find(genes.Chromosome == c & ...
			~any(isnan(genes.Position), 2));
		gene_map = false(1, organism.Chromosomes.Length(c));
		for g = chr_genes'
			gene_map(genes.Position(g, 1):genes.Position(g, 2)) = true;
		end
		intronic_reads(s) = intronic_reads(s) + sum(gene_map(offset(chr == c)));
	end
	
	intergenic_reads(s) = length(offset) - intronic_reads(s);
	
	unknown_reads(s) = total_reads(s) - sum(contaminant_reads(:, s)) - ...
		exonic_reads(s) - intronic_reads(s) - intergenic_reads(s);
end




% Write statistics into a spreadsheet.
fid = fopen([report_dir '/rnaseq_stats.txt'], 'W');
fprintf(fid, 'Sample\trRNA\tmtDNA\tExonic\tIntronic\tIntergenic\tUnknown\n');
for s = 1:S
	fprintf(fid, '%s\t%.1f%%\t%.1f%%\t%.1f%%\t%.1f%%\t%.1f%%\t%.1f%%\n', ...
		sample_ids{s}, ...
		sum(contaminant_reads(1:4, s)) / total_reads(s) * 100, ...
		contaminant_reads(5, s) / total_reads(s) * 100, ...
		exonic_reads(s) / total_reads(s) * 100, ...
		intronic_reads(s) / total_reads(s) * 100, ...
		intergenic_reads(s) / total_reads(s) * 100, ...
		unknown_reads(s) / total_reads(s) * 100);
end

fclose(fid);




