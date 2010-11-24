function [] = rename_dataset(name, new_name)

while new_name(end) == '/'
	new_name = new_name(1:end-1);
end

system(sprintf('mv %s/datasets/%s %s/datasets/%s', ppath, flatten_str(name), ...
	ppath, flatten_str(new_name)));

