
% This script is used for compiling all bacterial genomes into one FASTA
% file after downloading from the NCBI FTP the file
% ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.fna.tar.gz

% Author: Matti Annala <matti.annala@tut.fi>

function [] = compile_refseq_bacteria()

files = dir('.');
orig_dirs = files([files.isdir] == 1);
orig_dirs = {orig_dirs.name}';
orig_dirs = orig_dirs(3:end);

dirs = regexprep(orig_dirs, '_[A-Z0-9].+', '');
dirs = strrep(dirs, '_', ' ');

[species, idx] = unique(dirs, 'first');
orig_dirs = orig_dirs(idx);

fprintf(1, 'Found %d bacterial species...\n', length(species));

for k = 1:length(species)
	files = dir(orig_dirs{k});
	for f = 3:length(files)  % Skip '.' and '..'
		unix(sprintf('cat %s/%s >> bacteria.fa', ...
			orig_dirs{k}, files(f).name));
	end
end


