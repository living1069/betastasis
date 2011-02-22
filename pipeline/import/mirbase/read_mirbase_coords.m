function pre_mirnas = read_mirbase_coords(mirbase_gff, pre_mirnas)

fid = fopen(mirbase_gff);
data = textscan(fid, '%s %*s %*s %d %d %*s %s %*s %s', 'Delimiter', '\t');
fclose(fid);

chr = chromosome_sym2num(data{1});
left = data{2};
right = data{3};
strand = data{4};
extra = data{5};

pre_mirna_to_idx = containers.Map(pre_mirnas.Name, ...
	num2cell(1:length(pre_mirnas.Name)));

P = length(pre_mirnas.Name);
pre_mirnas.Chromosome = nan(P, 1);
pre_mirnas.Strand = repmat(' ', P, 1);
pre_mirnas.Position = nan(P, 2);

for k = 1:length(extra)
	tokens = regexpi(extra{k}, 'ID="(.+?)"', 'tokens');
	if length(tokens) ~= 1, continue, end
	
	token = tokens{1}; premir = token{1};
	if ~pre_mirna_to_idx.isKey(premir), continue, end
	
	idx = pre_mirna_to_idx(premir);
	
	pre_mirnas.Chromosome(idx, 1) = chr(k);
	pre_mirnas.Strand(idx, 1) = strand{k};
	pre_mirnas.Position(idx, 1:2) = [left(k) right(k)];
end

