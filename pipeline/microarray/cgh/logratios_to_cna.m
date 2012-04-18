
function cna = cgh_logratios_to_cna(logratios, test_meta, ref_meta, varargin)

global organism;

ploidy = cgh_ploidy(test_meta, ref_meta);

[P, S] = size(logratios);

if isfield(logratios, 'mean')
	for s = 1:S
		lr = logratios(:, s);
		
		for c = 1:length(organism.Chromosomes.Name)
			chr_segs = (results.chr == c);
			seg = struct;
			
			if ~isnan(ploidy(c, s))
				seg.Start = results.segstart(chr_segs);
				seg.End = results.segend(chr_segs);
			
				lr = results.segmean(chr_segs);
				seg.CNA = ploidy(c, s) * 2.^lr - ploidy(c, s);
				seg.CNA(seg.CNA < - 2) = -2;
			end
			
			segments.Chromosome{c, s} = seg;
		end
	end
	
	
elseif isfield(logratios, 'chromosome')
	for s = 1:S
		
	end
end

