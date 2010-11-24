function pooled = pool_rearrangements(rearrangements)

exon_pair_map = containers.Map();

pooled = struct('Exons', [], ...
                'ReadCount', []);

if isfield(rearrangements{1}, 'ReadSequences')
	pooled.ReadSequences = {};
end
				
for r = 1:length(rearrangements)
	sample = rearrangements{r};
	for k = 1:size(sample.Exons, 1)
		pair = sprintf('%d,%d', sample.Exons(k, 1), sample.Exons(k, 2));
		if ~exon_pair_map.isKey(pair)
			exon_pair_map(pair) = size(pooled.Exons, 1) + 1;
			pooled.Exons(end + 1, 1) = sample.Exons(k, 1);
			pooled.Exons(end, 2) = sample.Exons(k, 2);
			pooled.ReadCount(end + 1, 1) = 0;
		end
		idx = exon_pair_map(pair);
		pooled.ReadCount(idx) = pooled.ReadCount(idx) + sample.ReadCount(k);
		
		if isfield(sample, 'ReadSequences')
			pooled.ReadSequences{idx, pooled.ReadCount(idx)} = '';
			start = pooled.ReadCount(idx) - sample.ReadCount(k) + 1;
			pooled.ReadSequences(idx, start:pooled.ReadCount(idx)) = ...
				sample.ReadSequences(k, 1:sample.ReadCount(k));
		end
	end
end

