function [genes, affinities] = chipseq_genes(peak_file, refgene)

fprintf(1, 'Reading UCSC refGene annotations...\n');
refgene = read_mm9_refgene(refgene)

fprintf(1, 'Reading MACS peak file...\n');
fid = fopen(peak_file);
data = textscan(fid, '%s %d %d %*s %f', 'HeaderLines', 1, 'Delimiter', '\t');
fclose(fid);

chrs = unique(data{1});
chr_map = containers.Map(chrs, num2cell(1:length(chrs)));

refgene.Chromosome = cell2mat(chr_map.values(refgene.Chromosome));

peaks.Chromosome = cell2mat(chr_map.values(data{1}));
peaks.Start = data{2};
peaks.End = data{3};
peaks.Position = round((peaks.Start + peaks.End) / 2);
peaks.Intensity = data{4};

P = length(peaks.Start);

r = ~isnan(refgene.Chromosome) & ~isnan(refgene.Position(:, 1));
ranges.Gene = find(r);
ranges.Chromosome = refgene.Chromosome(r);
pos = refgene.Position(r, :);
strand = refgene.Strand(r);
plus = (strand == '+');
minus = (strand == '-');
ranges.Position(plus, :) = [pos(plus, 1)-10e3, pos(plus, 1)];
ranges.Position(minus, :) = [pos(minus, 2), pos(minus, 2)+10e3];

ranges

window_size = int32(10e6);

chr_ranges = repmat({{}}, 25, 1);
for c = 1:length(chr_ranges)
	chr_ranges{c} = cell(1, 100);
end

for r = 1:length(ranges.Gene)
	chr = ranges.Chromosome(r);
	wstart = idivide(ranges.Position(r, 1)-1, window_size) + 1;
	wend = idivide(ranges.Position(r, 2), window_size) + 1;
	
	for w = wstart:wend
		chr_ranges{chr}{w}(end+1) = r;
	end
end

affinities = zeros(length(refgene.Name), 1);

for p = 1:P
	chr = peaks.Chromosome(p);
	pos = peaks.Position(p);
	w = idivide(pos-1, window_size) + 1;
	
	nearby_ranges = chr_ranges{chr}{w};
	promoters = ranges.Position(nearby_ranges, :);
	
	overlapping = (promoters(:, 1) <= pos) & (promoters(:, 2) >= pos);
	overlapping = nearby_ranges(overlapping);
	
	affinities(ranges.Gene(overlapping)) = ...
		affinities(ranges.Gene(overlapping)) + peaks.Intensity(p);
end

genes = refgene.Name;









function refgene = read_mm9_refgene(filepath)

refgene = struct;

fid = fopen(filepath);
data = textscan(fid, '%s %s %s %s %d %d %d %d %d %s %s %s %s %s %s %s', ...
	'Delimiter', '\t');
fclose(fid);

chromosome = data{3};
strand = data{4};
left = data{5};
right = data{6};
gene_name = data{13};

refgene.Name = unique(gene_name);
N = length(refgene.Name);

refgene.Chromosome = cell(N, 1);
refgene.Position = nan(N, 2);
refgene.Strand = repmat(' ', N, 1);

invalid = false(N, 1);

gene_map = containers.Map(refgene.Name, num2cell(1:length(refgene.Name)));

gene_idx = cell2mat(gene_map.values(gene_name));

for k = 1:length(gene_idx)
	idx = gene_idx(k);
	
	if regexp(chromosome{k}, 'random')
		invalid(idx) = true;
	elseif isempty(refgene.Chromosome{idx})
		refgene.Chromosome{idx} = chromosome{k};
	elseif ~strcmp(refgene.Chromosome{idx}, chromosome{k})
		invalid(idx) = true;
	end
	
	if refgene.Strand(idx) == ' '
		refgene.Strand(idx) = strand{k};
	elseif ~strcmp(strand{k}, refgene.Strand(idx))
		invalid(idx) = true;
	end
	
	if left(k) > right(k)
		error 'Inverted coordinates.';
	end
	
	if isnan(refgene.Position(idx, 1))
		refgene.Position(idx, :) = [left(k) right(k)];
	elseif right(k) < refgene.Position(idx, 1) || ...
		left(k) > refgene.Position(idx, 2)
		invalid(idx) = true;
	else
		refgene.Position(idx, 1) = min(left(k), refgene.Position(idx, 1));
		refgene.Position(idx, 2) = max(right(k), refgene.Position(idx, 2));
	end
end

refgene.Name = refgene.Name(~invalid);
refgene.Strand = refgene.Strand(~invalid);
refgene.Chromosome = refgene.Chromosome(~invalid);
refgene.Position = refgene.Position(~invalid, :);

