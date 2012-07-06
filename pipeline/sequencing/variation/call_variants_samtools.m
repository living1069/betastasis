function [] = call_variants_samtools(alignments, out_file, varargin)

global organism;
chromosomes = organism.Chromosomes;

use_chr_prefix = true;
num_threads = 10;
remove_duplicates = false;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'chr.*prefix')
		use_chr_prefix = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'threads')
		num_threads = varargin{k+1}; continue;
	end
	
	if rx(varargin{k}, 'remove.*dup')
		remove_duplicates = varargin{k+1}; continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

S = length(alignments.url);

tmp = temporary('call_variants_samtools');

if use_chr_prefix
	ref_genome = '/data/csb/tools/bowtie2-indexes/homo_sapiens/hg19.fa';
else
	ref_genome = ['/data/csb/tools/bowtie2-indexes/homo_sapiens/' ...
		'hg19_no_chr_prefix.fa'];
end

variants = struct;
variants.meta = alignments.meta;

bam_file_paths = '';
for s = 1:length(alignments.url)
	bam_path = [alignments.url{s} '.bam'];
	if remove_duplicates
		bam_path = ['<(samtools rmdup -S ' bam_path ' -)'];
	end
	bam_file_paths = [bam_file_paths ' ' bam_path];
end

% -u = Output uncompressed BCF (fast for piping)
% -B = Disable BAQ (increases sensitivity)
% -D = Output per-sample read depths (DP)
mpileup_flags = '-uBD';

bcf_view_flags = '-vcg -p 0.1';

num_parallel = 10;
chr_status = zeros(1, length(chromosomes.Name));
while ~all(chr_status == 2)
	for c = [23 1:22 24]   % X chromosome is large, so we do it early
		chr_out_prefix = [tmp chromosomes.Name{c}];
		
		if chr_status(c) == 0 && sum(chr_status == 1) < num_parallel
			
			if use_chr_prefix
				chr_symbol = ['chr' chromosomes.Name{c}];
			else
				chr_symbol = chromosomes.Name{c};
			end
			
			daemonize = ' &';
			if num_threads == 1, daemonize = ''; end			
			
			unix(sprintf(['samtools mpileup %s -f %s -r %s %s | ' ...
				'bcftools view %s - | gzip -c > %s.incomplete.vcf.gz && ' ...
				'mv %s.incomplete.vcf.gz %s.vcf.gz %s'], ...
				mpileup_flags, ref_genome, chr_symbol, bam_file_paths, ...
				bcf_view_flags, chr_out_prefix, chr_out_prefix, ...
				chr_out_prefix, daemonize));
			chr_status(c) = 1;
		elseif chr_status(c) == 1 && exist([chr_out_prefix '.vcf.gz'])
			chr_status(c) = 2;
		end
	end
	pause(10);
end

chr_str = sprintf([' ' tmp '%s.vcf.gz'], chromosomes.Name{:});
unix(sprintf('cat %s > %s.gz', chr_str, out_file));

