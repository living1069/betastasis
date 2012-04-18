
% Author: Matti Annala <matti.annala@tut.fi>

function reads = sra_to_fasta(root, varargin)

delete_originals = true;

sra_files = find_files('.*\.sra');

for f = 1:length(sra_files)
	out_dir = regexprep(sra_files{f}, '/[^/]*', '');
	
	[status, out] = unix(sprintf(['%s/tools/sratoolkit.2.1.0-centos_linux64/'...
		'fastq-dump --fasta -O %s %s'], ppath, out_dir, sra_files{f}));
	if status ~= 0, error('SRA to FASTQ dump failed:\n%s\n', out); end
	
	if delete_originals, delete(sra_files{f}); end
end


