function [] = call_variants_samtools(alignments, out_file)

global organism;
chromosomes = organism.Chromosomes;

S = length(alignments.url);

tmp = temporary('call_variants_samtools');

variants = struct;
variants.meta = alignments.meta;

bam_file_paths = '';
for s = 1:length(alignments.url)
	bam_file_paths = [bam_file_paths ' ' alignments.url{s} '.bam'];
end

% -u = Output uncompressed BCF (fast for piping)
% -B = Disable BAQ (increases sensitivity)
% -D = Output per-sample read depths (DP)
mpileup_flags = '-uBD -f ~/tools/bowtie2-indexes/homo_sapiens/2009/genome.fa';

bcf_view_flags = '-vcg -p 0.1';

num_parallel = 8;
chr_status = zeros(1, length(chromosomes.Name));
while ~all(chr_status == 2)
	for c = 1:length(chr_status)
		if chr_status(c) == 0 && sum(chr_status == 1) < num_parallel
			chr_out_prefix = [tmp chromosomes.Name{c}];
			unix(sprintf(['samtools mpileup %s -r chr%s %s | ' ...
				'bcftools view %s - | gzip -c > %s.incomplete.vcf.gz && ' ...
				'mv %s.incomplete.vcf.gz %s.vcf.gz'], ...
				mpileup_flags, chromosomes.Name{c}, bam_file_paths, ...
				bcf_view_flags, chr_out_prefix, chr_out_prefix, ...
				chr_out_prefix));
			chr_status(c) = 1;
		elseif chr_status(c) == 1 && exist([tmp chromosomes.Name{c} '.vcf.gz'])
			chr_status(c) = 2;
		end
	end
	pause(10);
end

chr_str = '';
for c = 1:length(chromosomes.Name)
	chr_str = [chr_str ' ' chromosomes.Name{c} '.vcf.gz'];
end

unix(sprintf('cat %s > %s.gz', chr_str, out_file));


