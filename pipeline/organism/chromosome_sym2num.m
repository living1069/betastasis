function chrnum = chromosome_sym2num(symbol)

if iscell(symbol)
	for k = 1:numel(symbol)
		chrnum(k) = sym2num(symbol{k});
	end
elseif ischar(symbol)
	chrnum = sym2num(symbol);
end

return;





function chrnum = sym2num(symbol)

global organism;

symbol = strrep(symbol, ' ', '');

idx = find(strcmp(symbol, organism.Chromosomes.Name));
if length(idx) > 0
	chrnum = idx(1);
	return;
end

if length(symbol) > 3 && strcmp(symbol(1:3), 'chr')
	idx = find(strcmp(symbol(4:end), organism.Chromosomes.Name));
	if length(idx) > 0
		chrnum = idx(1);
		return;
	end
end

if strcmpi(symbol, 'Unknown') || isempty(symbol)
	chrnum = 0;
	return;
end

chrnum = NaN;
return;
