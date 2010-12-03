function [] = select_organism(name, version)

name(1) = upper(name(1));

orgpath = [ppath '/organisms/' flatten_str(name) '/' flatten_str(version)];
if exist(orgpath) ~= 7
	fprintf(1, 'Could not find the specified organism.\n');
	return;
end

global organism;
organism = Organism(name, version);

