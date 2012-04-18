function [] = dnaseq_stats(reads, varargin)

global organism;

%load ~/organisms/homo_sapiens/ucsc_hg19/refgene.mat;
%genes = refgene;
%transcripts = refgene.Transcripts;

load ~/organisms/homo_sapiens/ensembl_59/ensembl.mat;
genes = ensembl.Genes;
transcripts = ensembl.Transcripts;

seq_files = seq_resource_files(reads);

chr_ranges = zeros(length(organism.Chromosomes.Name), 2);
for k = 1:length(organism.Chromosomes.Name)
	chr_ranges(k, :) = [1 + sum(organism.Chromosomes.Length(1:k-1), 1), ...
		sum(organism.Chromosomes.Length(1:k), 1)];
end

gene_map = false(1, chr_ranges(end, 2));
exon_map = false(1, chr_ranges(end, 2));
for ts = 1:length(transcripts.Name)
	if isnan(transcripts.Chromosome(ts)), continue, end
	
	exons = transcripts.Exons{ts};
	range = [exons(1, 1) exons(end, 2)] + ...
		chr_ranges(transcripts.Chromosome(ts), 1) - 1;
	if range(1) > range(2), range = range(end:-1:1); end
	gene_map(range(1):range(2)) = true;
	
	for e = 1:size(exons, 1)
		range = exons(e, :) + chr_ranges(transcripts.Chromosome(ts), 1) - 1;
		if range(1) > range(2), range = range(end:-1:1); end
		exon_map(range(1):range(2)) = true;
	end
end

run_ends = [ find(gene_map(1:end-1) ~= gene_map(2:end)), length(gene_map) ];
run_lengths = diff([0, run_ends]);

gene_pos = [];
for r = 1:length(run_ends)
	if gene_map(run_ends(r)) == 0, continue, end
	gene_pos(end+1, :) = [run_ends(r-1)+1, run_ends(r)];
end

upstream_2kb_map = false(1, chr_ranges(end, 2));
for g = 1:size(gene_pos, 1)
	range = [gene_pos(g, 1)-2000, gene_pos(g, 1)-1];
	upstream_2kb_map(range(1):range(2)) = true;
end
upstream_2kb_map = upstream_2kb_map & (~gene_map);

fprintf(1, 'Gene coverage of genome: %.1f%%\n', ...
	100 * sum(gene_map) / length(gene_map));
fprintf(1, 'Exon coverage of genome: %.1f%%\n', ...
	100 * sum(exon_map) / length(exon_map));

alignments = cell(length(seq_files), 1);

%load ~/al_tmp.mat;
%alignments = al;

for seq_file = 1:length(seq_files)
	alignments{seq_file} = align_reads(seq_files{seq_file}, 'genome', ...
		varargin{:}, 'AllowAlignments', 1, 'Columns', 'target,offset,sequence');
	al = alignments{seq_file};

	chromosomes = zeros(length(al.Target), 1);
	for k = 1:length(al.Target)
		chromosomes(k) = chromosome_sym2num(al.Target{k});
	end
	
	seqlens = zeros(length(al.Sequence), 1);
	for k = 1:length(seqlens)
		seqlens(k) = length(al.Sequence{k});
	end
	
	offsets = double(al.Offset);
	offsets(1:10)
	offsets = round((offsets + seqlens) / 2);
	offsets = offsets + chr_ranges(chromosomes, 1) - 1;
	offsets(1:10)
	chromosomes(1:10)
	chr_ranges(1:10, 1)
	
	N = length(offsets);
	gene_hits = 0;
	exon_hits = 0;
	upstream_2kb_hits = 0;
	
	for k = 1:N
		if gene_map(offsets(k)), gene_hits = gene_hits + 1; end
		if exon_map(offsets(k)), exon_hits = exon_hits + 1; end
		if upstream_2kb_map(offsets(k))
			upstream_2kb_hits = upstream_2kb_hits + 1;
		end
	end
	
	intron_hits = gene_hits - exon_hits;
	
	fprintf(1, 'Read distribution by genomic feature:\n');
	fprintf(1, '- reads from exons: %d (%.1f%%)\n', exon_hits, ...
		exon_hits / N * 100);
	fprintf(1, '- reads from introns: %d (%.1f%%)\n', intron_hits, ...
		intron_hits / N * 100);
	fprintf(1, '- reads from gene upstream regions (2kb): %d (%.1f%%)\n', ...
		upstream_2kb_hits, upstream_2kb_hits / N * 100);
end


