function darned = read_darned(filepath)
	
global organism;

darned = struct;

fid = fopen(filepath);
data = textscan(fid, '%s%d%s%s%s%s%s%s%s%s%s%s%s', 'HeaderLines', 1, ...
	'Delimiter', '\t');
fclose(fid);

N = length(data{1});

darned.Chromosome = chromosome_sym2num(data{1})';
darned.Offset = data{2};

darned.Strand = repmat(' ', N, 1);
darned.RefBase = repmat(' ', N, 1);
darned.EditBase = repmat(' ', N, 1);

for k = 1:N
	darned.Strand(k) = data{3}{k};
	darned.RefBase(k) = data{4}{k};
	darned.EditBase(k) = data{5}{k};
end

darned.Tissue = data{10};

