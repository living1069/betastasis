function m = rx(x, regex)

m = false(length(x), 1);
for k = 1:length(x)
	if regexpi(x{k}, regex), m(k) = true; end
end

