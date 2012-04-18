function [path, dir] = path_strip_dir(path)
pos = find(path == '/');
dir = '.';
if isempty(pos), return, end
dir = path(1:pos(end));
path = path(pos(end)+1:end);

