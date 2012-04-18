function pooled = pool_rearrangements(rearrangements)

exon_pair_map = containers.Map;

% Calculate an upper limit for the number of pooled fusion candidates.
prealloc = 0;
for s = 1:length(rearrangements)
	prealloc = prealloc + size(rearrangements{s}.Exons, 1);
end

pooled = struct;
pooled.Exons = nan(prealloc, 2);
pooled.ReadCount = zeros(prealloc, 1);

F = 0;
progress = Progress;
				
for r = 1:length(rearrangements)
	sample = rearrangements{r};
	
	for k = 1:size(sample.Exons, 1)
		pair = sprintf('%d,%d', sample.Exons(k, 1), sample.Exons(k, 2));
		if ~exon_pair_map.isKey(pair)
			F = F + 1;
			exon_pair_map(pair) = F;
			pooled.Exons(F, :) = sample.Exons(k, :);
		end
		idx = exon_pair_map(pair);
		pooled.ReadCount(idx) = pooled.ReadCount(idx) + sample.ReadCount(k);
	end
	
	progress.update(r / length(rearrangements));
end

pooled.Exons = pooled.Exons(1:F, :);
pooled.ReadCount = pooled.ReadCount(1:F);

