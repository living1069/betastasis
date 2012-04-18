function [] = fprintf_cell(fid, x, delim)

if ~isempty(x)
	fprintf(fid, sprintf('%s%s', '%s', delim), x{:});
	%for s = 1:length(x)-1, fprintf(fid, '%s%s', x{s}, delim); end
	fprintf(fid, '%s', x{end});
end

