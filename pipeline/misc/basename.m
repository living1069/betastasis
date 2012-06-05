function name = basename(path)
pos = find(path == '/');
if isempty(pos)
	name = path; return;
end
name = path(pos(end)+1:end);

