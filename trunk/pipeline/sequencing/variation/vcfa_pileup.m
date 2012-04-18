
function [] = vcfa_pileup(vcfa_file, variant_indices, varargin)

global organism;
chromosomes = organism.Chromosomes;

tmp = temporary('vcfa_pileup');

bam_root = '.';

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'bam.*root')
		bam_root = varargin{k+1}; continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

files = find_files('.bam$', bam_root);
files = files(~rx(files, '.vcfa.bam'));

% Read variant coordinates from the VCFA file and calculate flanking regions.
data = readtable(vcfa_file);
chr = data{1}(variant_indices);
offset = str2double(data{2}(variant_indices));

for k = 1:length(chr)
	regions{k} = sprintf('%s:%d-%d', chr{k}, offset(k)-150, offset(k)+150);
end

for s = 1:length(files)
	filtered_bams{s} = sprintf('%s/%d.bam', tmp, s);
end

for s = 1:length(files)
	[status, out] = unix(sprintf('samtools view -b %s %s > %s', files{s}, ...
		sprintf('%s ', regions{:}), filtered_bams{s}));
	if status ~= 0, error('samtools view failed:\n%s\n', out); end
end

[status, out] = unix(sprintf('samtools merge -r %s %s', ...
	[vcfa_file, '.bam'], sprintf('%s ', filtered_bams{:})));

	

