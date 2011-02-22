function region_cna = feature_cna(segments, cgh_probesets, regions, varargin)

cna = cn_seg_expand(segments, cgh_probesets);

S = size(samples, 2);
R = length(regions.Chromosome);
region_cna = zeros(R, S);

for r = 1:R
	ps = find(cgh_probesets.Chromosome == regions.Chromosome(r) & ...
		cgh_probesets.Offset >= regions.Position(1) & ...
		cgh_probesets.Offset <= regions.Position(2));
	
	region_cna(r, :) = nanmedian(cna(ps, :), 1);
end

return;




for f = 1:length(features.Spec)
	tokens = regexpi(features.Spec{f}, ...
		'chr(.+?):\s*(\d+)-(\d+)-(\d+)-(\d+)', 'tokens');
	if length(tokens) == 1
		token = tokens{1};
		chr = chromosome_sym2num(token{1});
		if str2double(token{2}) >= str2double(token{3}) || ...
			str2double(token{3}) >= str2double(token{4}) || ...
			str2double(token{4}) >= str2double(token{5})
			error 'Invalid CNA feature region.';
		end
		
		mid_idx = find(cgh_probesets.Chromosome == chr & ...
			cgh_probesets.Offset >= str2double(token{3}) & ...
			cgh_probesets.Offset <= str2double(token{4}));
		outer_idx = find(cgh_probesets.Chromosome == chr & ...
			cgh_probesets.Offset >= str2double(token{2}) & ...
			cgh_probesets.Offset <= str2double(token{5}));
		outer_idx = setdiff(outer_idx, mid_idx);
			
		feature_cna(f, :) = nanmedian(cna(mid_idx, :), 1) - ...
			median(cna(outer_idx, :), 1);
		continue;
	end
	
	tokens = regexpi(features.Spec{f}, ...
		'chr(.+?):\s*(\d+)-(\d+)', 'tokens');
	if length(tokens) == 1
		token = tokens{1};
		chr = chromosome_sym2num(token{1});
		if str2double(token{2}) >= str2double(token{3})
			error 'Invalid CNA feature region.';
		end
		
		idx = find(cgh_probesets.Chromosome == chr & ...
			cgh_probesets.Offset >= str2double(token{2}) & ...
			cgh_probesets.Offset <= str2double(token{3}));
			
		feature_cna(f, :) = nanmedian(cna(idx, :), 1);
	end
end


