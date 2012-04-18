
function [] = store(name, input_data)

global pipeline_config;
if isempty(pipeline_config.Repositories)
	error 'No repositories have been configured.';
end

root = pipeline_config.Repositories{1}.URL;   % FIXME

check_dataset_path(name);

name = flatten_str(name);
path = [root '/' name];
[~, ~, ~] = mkdir(path);

metadata = input_data.Meta;
input_data = rmfield(input_data, 'Meta');

S = 0;
fields = fieldnames(input_data);
for k = 1:length(fields)
	f = getfield(input_data, fields{k});
	if (isnumeric(f) || iscell(f)) && numel(f) > 0
		S = size(f, 2);
		break;
	end
end

if S == 0, error 'Could not determine data set dimensions.'; end
	
fprintf(1, 'Data set contains %d samples.\n', S);

for k = 1:S
	s = struct;
	for n = 1:length(fields)
		f = getfield(input_data, fields{n});
		if size(f, 2) == S
			eval(sprintf('s.%s = f(:, %d);', fields{n}, k));
		else
			error('create_dataset() could not handle field "%s".\n', ...
				fields{n});
		end
	end
	metadata.Resource{k, 1} = create_sample(name, s);
end

save([path '/metadata.mat'], 'metadata');

