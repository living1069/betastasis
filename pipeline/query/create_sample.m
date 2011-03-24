function resource_url = create_sample(dataset, data)

global pipeline_config;
root = pipeline_config.Repositories{1}.URL;   % FIXME

if exist([root '/' flatten_str(dataset)]) ~= 7
	[~, ~] = mkdir([root '/' flatten_str(dataset)]);
end

if isfield(data, 'Raw') && length(data.Raw) == 1
	resource_url = '';
	for f = 1:length(data.Raw{1}.Paths)
		path = data.Raw{1}.Paths{f};
		resource = gen_uuid();
		if strcmpi(path(end-2:end), '.gz')
			resource = [resource '.gz'];
		end
		copyfile(data.Raw{1}.Paths{f}, ...
			[root '/' flatten_str(dataset) '/' resource]);
		resource_url = [resource_url resource ';'];
	end
	resource_url = resource_url(1:end-1);  % Remove the last semicolon.
else
	resource_url = gen_uuid();
	save([root '/' flatten_str(dataset) '/' resource_url], 'data');
end

