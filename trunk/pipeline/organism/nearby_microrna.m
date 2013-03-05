function features = nearby_microrna(features)

max_distance = 100e3;

data = readtable( ...
	'~/organisms/homo_sapiens/mirbase_18/miRNA_26_3_2012_parsed.txt', ...
	'Header', false, 'Numeric', 2:3);
chr = chromosome_sym2num(data{1});
pos = [data{2} data{3}];
strand = char(data{4});
mirna = data{6};
tss_pos = nan(size(chr));
tss_pos(strand == '+') = pos(strand == '+', 1);
tss_pos(strand == '-') = pos(strand == '-', 1);

F = length(features.chromosome);

for k = 1:F
	dist = features.position(k) - tss_pos;
	dist = dist .* ((strand == '+') - 0.5) * 2;
	nearby = find(chr == features.chromosome(k) & abs(dist) < max_distance);

	features.nearby_mirna{k, 1} = mirna(nearby);
	features.nearby_mirna_dist{k, 1} = dist(nearby);
	
	% Discard identical entries.
	tmp = cell(length(features.nearby_mirna{k}), 1);
	for n = 1:length(tmp)
		tmp{n} = sprintf('%s (%d)', features.nearby_mirna{k}{n}, ...
			features.nearby_mirna_dist{k}(n));
	end
	[~, uniq] = unique(tmp);
	features.nearby_mirna{k} = features.nearby_mirna{k}(uniq);
	features.nearby_mirna_dist{k} = features.nearby_mirna_dist{k}(uniq);
	
	% Sort according to distance.
	[~, order] = sort(abs(features.nearby_mirna_dist{k}));
	features.nearby_mirna{k} = features.nearby_mirna{k}(order);
	features.nearby_mirna_dist{k} = features.nearby_mirna_dist{k}(order);
end


