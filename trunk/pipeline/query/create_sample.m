function resource_url = create_sample(dataset, data)

global pipeline_config;
root = pipeline_config.Repositories{1}.URL;   % FIXME

if exist([root '/' flatten_str(dataset)]) ~= 7
	system(['mkdir -p ' root '/' flatten_str(dataset)]);
end

resource_url = gen_uuid();
save([root '/' flatten_str(dataset) '/' resource_url], 'data');

