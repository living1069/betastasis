function [] = call_variants_varscan(alignments, out_file)

global organism;
chromosomes = organism.Chromosomes;

tmp = temporary('call_variants_varscan');

min_variant_allele_freq = 0.1;

S = length(alignments.url);

variants = struct;
variants.meta = alignments.meta;

bam_file_paths = '';
for s = 1:length(alignments.url)
	bam_file_paths = [bam_file_paths ' ' alignments.url{s} '.bam'];
end

% -u = Output uncompressed BCF (fast for piping)
% -B = Disable BAQ (increases sensitivity)
% -D = Output per-sample read depths (DP)
mpileup_flags = sprintf( ...
	'-B -f %s/tools/bowtie/indexes/homo_sapiens/2011/genome.fa', ppath);

varscan_flags = ['--min-coverage 10 --min-reads2 5 --min-var-freq 0.10 ' ...
	'--min-freq-for-hom 0.80 --variants'];
	
% Start multiple VarScan runs, each looking at a different chromosome.
for c = 1:length(chromosomes.Name)
	chr_out_prefix = [tmp chromosomes.Name{c}];
	[status, out] = unix(sprintf(['samtools mpileup %s -r chr%s %s | ' ...
		'java -jar /home/csbgroup/tools/VarScan.v2.2.10.jar mpileup2cns ' ...
		'--output-vcf %s > %s.incomplete.vcf && ' ...
		'mv %s.incomplete.vcf %s.vcf &'], ...
		mpileup_flags, chromosomes.Name{c}, bam_file_paths, ...
		varscan_flags, chr_out_prefix, chr_out_prefix, chr_out_prefix));
	if status ~= 0, error('samtools mpileup | VarScan failed:\n%s', out); end
end

for c = 1:length(chromosomes.Name)
	while ~exist([tmp chromosomes.Name{c} '.vcf'])
		pause(10);
	end
end

chr_str = '';
for c = 1:length(chromosomes.Name)
	chr_str = [chr_str ' chr' chromosomes.Name{c}];
end

unix(sprintf('cat %s > %s', chr_str, out_file));


