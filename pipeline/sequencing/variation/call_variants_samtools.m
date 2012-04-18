function [] = call_variants_samtools(alignments, out_file)

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
	'-uBD -f %s/tools/bowtie/indexes/homo_sapiens/2011/genome.fa', ppath);

bcf_view_flags = '-vcg -p 0.01';

vcf_filter_flags = sprintf(' -d %d', 5*S);   % At least 5 reads per sample


%if groups is not None:
%	alignments = alignments[hcat(groups[0], groups[1])]
%	system('%s/tools/samtools/samtools mpileup %s %s | '
%		'%s/tools/samtools/bcftools/bcftools view %s -1 %d - | '
%		'%s/tools/samtools/bcftools/vcfutils.pl varFilter %s > %s ' %
%		(pypette.path, mpileup_flags, bam_file_paths,
%		pypette.path, bcf_view_flags, len(group1),
%		pypette.path, vcf_filter_flags, out_file))

[status, out] = unix(sprintf(['%s/tools/samtools/samtools mpileup %s %s | ' ...
	'%s/tools/samtools/bcftools/bcftools view %s - | ' ...
	'%s/tools/samtools/bcftools/vcfutils.pl varFilter %s > %s '], ...
	ppath, mpileup_flags, bam_file_paths, ppath, bcf_view_flags, ...
	ppath, vcf_filter_flags, out_file));
if status ~= 0, error('samtools mpileup | bcftools failed:\n%s', out); end


