function [] = list_organisms()

fprintf(1, 'List of available organisms:\n');

files = dir([ppath '/organisms']);
for k = 1:length(files)
	if ~files(k).isdir, continue, end
	
	name = files(k).name;
	if strcmp(name, '..') || strcmp(name, '.'), continue, end
	
	name(1) = upper(name(1));
	name = strrep(name, '_', ' ');
	fprintf(1, '- %s\n', name);
	
	%vfiles = dir([ppath '/organisms/' files(k).name]);
	%for n = 1:length(vfiles)
	%	vname = vfiles(n).name;
	%	if isempty(regexp(vname, '\.mat$')), continue, end
	%	vname = vname(1:end-4);
	%	vname = strrep(vname, '_', ' ');
	%	fprintf(1, '  * %s\n', vname);
	%end
end

