function [] = call_variants_varscan(alignments, out_file)

global organism;
chromosomes = organism.Chromosomes;

tmp = temporary('call_variants_varscan');
%tmp = regexprep(out_file, '[^/]*$', '');

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
mpileup_flags = '-B -f ~/tools/bowtie2-indexes/homo_sapiens/2009/genome.fa';

varscan_flags = ['--min-coverage 10 --min-reads2 5 --min-var-freq 0.20 ' ...
	'--min-freq-for-hom 0.80 --variants'];
	
% Start multiple VarScan runs, each looking at a different chromosome.
num_parallel = 8;
chr_status = zeros(1, length(chromosomes.Name));
chr_status(:) = 2; chr_status(23) = 0;  % FIXME!!
while ~all(chr_status == 2)
	for c = 1:length(chr_status)
		if chr_status(c) == 0 && sum(chr_status == 1) < num_parallel
			chr_out_prefix = [tmp chromosomes.Name{c}];
			unix(sprintf([ ...
				'samtools mpileup %s -r chr%s %s | ' ...
				'java -jar ~/tools/VarScan.v2.2.11.jar mpileup2cns ' ...
				'--output-vcf %s > %s.incomplete.vcf && ' ...
				'mv %s.incomplete.vcf %s.vcf &'], ...
				mpileup_flags, chromosomes.Name{c}, bam_file_paths, ...
				varscan_flags, chr_out_prefix, chr_out_prefix, chr_out_prefix));
			chr_status(c) = 1;
		elseif chr_status(c) == 1 && exist([tmp chromosomes.Name{c} '.vcf'])
			chr_status(c) = 2;
		end
	end
	pause(10);
end

chr_str = '';
for c = 1:length(chromosomes.Name)
	chr_str = [chr_str ' chr' chromosomes.Name{c}];
end

unix(sprintf('cat %s > %s', chr_str, out_file));


