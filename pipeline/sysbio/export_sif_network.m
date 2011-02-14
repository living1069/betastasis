
% Author: Matti Annala <matti.annala@tut.fi>

function [] = export_sif_network(links, sif_file, varargin)

global organism;
genes = organism.Genes;

for k = 1:2:length(varargin)
	error('Unrecognized option "%s".', varargin{k});
end

tfs = unique(links.Genes(:, 1));

fid = fopen(sif_file, 'W');
for k = 1:length(tfs)
	tf = tfs(k);
	fprintf(fid, '%s mi', genes.Name{tf});
	targets = links.Genes(links.Genes(:, 1) == tf, 2);
	for t = 1:length(targets)
		fprintf(fid, ' %s', genes.Name{targets(t)});
	end
	fprintf(fid, '\n');
end
fclose(fid);

