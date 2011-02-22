
% Author: Matti Annala <matti.annala@tut.fi>

function [] = export_network_xls(links, xls_file, varargin)

global organism;
genes = organism.Genes;

feature_names = [];
subnetwork = true(max(max(links.Genes)), 1);

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'FeatureNames')
		feature_names = varargin{k+1};
		if ~iscellstr(feature_names)
			error 'Feature names must be provided as a cell array of strings.';
		end
		continue;
	end
	
	if strcmpi(varargin{k}, 'Subnetwork')
		subnetwork = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

if ~isempty(feature_names)
	if max(max(links.Genes)) > length(feature_names)
		error 'Provided feature names do not match with network nodes.';
	end
else
	feature_names = genes.Name;
	if max(max(links.Genes)) > length(feature_names)
		error(['No feature names were provided, and the network ' ...
			'has more nodes than there are genes for the organism.']);
	end
end

fid = fopen(xls_file, 'W');

tfs = unique(links.Genes(:, 1));

	
	
	
for k = 1:length(tfs)
	if subnetwork(tfs(k)) == false, continue, end

	tf = tfs(k);
	fprintf(fid, '%s', feature_names{tf});
	targets = links.Genes(links.Genes(:, 1) == tf, 2);
	for t = 1:length(targets)
		fprintf(fid, '\t%s', feature_names{targets(t)});
	end
	fprintf(fid, '\n');
end
fclose(fid);

