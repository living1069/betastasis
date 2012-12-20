function cna = cn_seg_expand(segs, probesets)

global organism;

S = size(segs.chromosome, 2);
cna = zeros(length(probesets.Offset), S);

chr_probes = nan(length(organism.Chromosomes.Name), 2);

for chr = 1:length(organism.Chromosomes.Name)
	idx = find(probesets.Chromosome == chr);
	if ~isempty(idx)
		chr_probes(chr, :) = [min(idx) max(idx)];
	end
end

for s = 1:S
	for chr = 1:size(segs.chromosome, 1)
		if any(isnan(chr_probes(chr, :))), continue, end
			
		probe_offsets = probesets.Offset(chr_probes(chr, 1):chr_probes(chr, 2));
		chr_segs = segs.chromosome{chr, s};
		if isempty(chr_segs) || ~isfield(chr_segs, 'Start'), continue, end
		
		for seg = 1:length(chr_segs.Start)
			seg_probes = chr_probes(chr, 1) - 1 + ...
				find(probe_offsets >= chr_segs.Start(seg) & ...
				probe_offsets <= chr_segs.End(seg));
			cna(seg_probes, s) = chr_segs.CNA(seg);
		end
	end
end

