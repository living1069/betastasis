function path = path_strip_dir(path)
pos = find(path == '/');
if isempty(pos), return, end
path = path(pos(end)+1:end);

