
% CREATE_DATASET    Store data as a named dataset in a pipeline repository
%
%    CREATE_DATASET(NAME, DATA) creates a named dataset in the current
%    default repository. The dataset contains the data stored in argument DATA,
%    which must be structured according to the pipeline data model
%    (i.e. must have metadata under field "Meta").
%
%    The dataset name NAME can be prefixed with directory names separated by
%    slash ('/') characters. This allows datasets to be stored in a hierarchical
%    tree within repositories. The NAME can also be prefixed with a repository
%    identifier followed by a colon. This allows datasets to be stored in
%    repositories other than the current default repository. Here is an example:
%
%        create_dataset('remote repo:cancer/glioma/mutations', mutations)
%
%    See also REMOVE_DATASET, QUERY.

function [] = create_dataset(name, input_data)

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

N = 0;
fields = fieldnames(input_data);
for k = 1:length(fields)
	f = getfield(input_data, fields{k});
	if (isnumeric(f) || iscell(f)) && numel(f) > 0
		N = size(f, 2);
		break;
	end
end

if N == 0, error 'Could not determine data set dimensions.'; end
	
fprintf(1, 'Data set appears to consist of %d samples.\n', N);

for k = 1:N
	s = struct;
	for n = 1:length(fields)
		f = getfield(input_data, fields{n});
		if size(f, 2) == N
			eval(sprintf('s.%s = f(:, %d);', fields{n}, k));
		else
			error('create_dataset() could not handle field "%s".\n', ...
				fields{n});
		end
	end
	metadata.Resource{k, 1} = create_sample(name, s);
end

save([path '/metadata.mat'], 'metadata');

