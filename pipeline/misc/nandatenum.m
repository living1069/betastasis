function num = nandatenum(str)
N = numel(str);
num = nan(size(str));
empty = cellfun(@isempty, str);
num(~empty) = datenum(str(~empty));

