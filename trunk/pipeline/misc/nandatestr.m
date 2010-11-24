function str = nandatestr(n)
N = numel(n);
str = repmat({'-'}, N, 1);
nnan_idx = find(~isnan(n(:)));
dates = datestr(n(nnan_idx));
for k = 1:length(nnan_idx)
	str{nnan_idx(k)} = dates(k, :);
end

