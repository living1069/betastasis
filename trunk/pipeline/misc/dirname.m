function name = dirname(path)
pos = find(path == '/');
if isempty(pos)
	name = '.';
	return;
end
name = path(1:pos(end)-1);

