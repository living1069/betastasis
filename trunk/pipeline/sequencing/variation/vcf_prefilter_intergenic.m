
function [] = vcf_prefilter_intergenic(vcf_file)

tmp = temporary('vcf_prefilter_intergenic');
out_prefix = [tmp 'variants'];

% Command used to convert the VCF file into an ANNOVAR compatible format.
% Note that the 'sed' command is used to fix a bug in convert2annovar.pl
convert_cmd = sprintf(['convert2annovar.pl --format vcf4 --includeinfo ' ...
	'--allallele <(gunzip -c %s) | sed "s/\t\t/\t/"'], vcf_file);

humandb_dir = '~/tools/annovar*/humandb';
	
unix(sprintf(['annotate_variation.pl --buildver hg19 --geneanno ' ...
	'--outfile %s <(%s) %s'], out_prefix, convert_cmd, humandb_dir));

unix(sprintf(['grep -v intergenic %s.variant_function > ' ...
	'%s_intragenic.variant_function'], out_prefix, out_prefix));
	
unix(sprintf('cut -f 8- %s_intragenic.variant_function | gzip -c > %s', ...
	out_prefix, strrep(vcf_file, '.vcf.gz', '.intragenic.vcf.gz')));

