function cosmic = import_cosmic(file)

[data, headers] = readtable(file);

cosmic = struct;

grch37_pos = data{rx(headers, 'grch37')};
pubmed = data{rx(headers, 'pubmed')};
cosmic.map = containers.Map( ...
	regexprep(grch37_pos, '^(.+:\d+)-\d+$', 'chr$1'), pubmed);


