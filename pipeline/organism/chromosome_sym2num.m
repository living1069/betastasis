function chrnum = chromosome_sym2num(symbol)

global organism;

if ischar(symbol), symbol = { symbol }; end

symbol = strrep(symbol, 'chr', '');
symbol = strrep(symbol, ' ', '');

tmp = strcmp('MT', symbol);
symbol(tmp) = repmat({'M'}, sum(tmp), 1);

[~, chrnum] = ismember(symbol, organism.Chromosomes.Name);
chrnum(chrnum == 0) = NaN;

