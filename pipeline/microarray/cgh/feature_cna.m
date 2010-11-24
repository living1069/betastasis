function feature_cna = feature_cna(samples, refs, cgh_probesets, features, ...
	varargin)

cna = cna_from_cgh(samples, refs, cgh_probesets, varargin{:});

S = size(samples, 2);

feature_cna = zeros(length(features.Name), S);

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


