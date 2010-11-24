function data = realize(qset)

global pipeline_config;

data = struct;

if ~isfield(qset, 'Resource')
	error 'realize() only works on query sets.';
end

N = length(qset.Resource);
if N == 0, error('Query data set is empty.\n'); end

fprintf(1, 'Progress: 00%%');
progress = 0;

resources = strcat(qset.Resource, '.mat');
resource_files = resources;

% Find a handle to the repository where the resources are stored and
% cache (download) all referred data resources.
repo_ids = {};
for k = 1:length(resources)
	colon = find(resources{k} == ':');
	if length(colon) ~= 1, error 'Invalid resource identifier.'; end
	repo_ids{end + 1, 1} = resources{k}(1:colon-1);
	resource_files{k} = resources{k}(colon+1:end);
end

repo_name = unique(repo_ids);
if length(repo_name) > 1
	error 'Realizing multi-repo datasets is not yet supported.';
end
repo_name = repo_name{1};

for r = 1:length(pipeline_config.Repositories)
	repo = pipeline_config.Repositories{r};
	if strcmp(repo_name, repo.Name), break, end
	if r == length(pipeline_config.Repositories)
		error('Repository "%s" was not found in configuration.', repo_name);
	end
end

repo.cache(resource_files);

% Load the first data resource and extend its column vectors into matrices
% in preparation for loading the other data resources.
s = fetch_resource(resource_files{1}, repo);
fields = fieldnames(s);

for k = 1:length(fields)
	f = getfield(s, fields{k});
	if isnumeric(f) && size(f, 2) == 1
		eval(['data.' fields{k} ' = f; data.' fields{k} '(end, N) = 0;']);
	elseif iscell(f) && size(f, 2) == 1
		data = setfield(data, fields{k}, cat(2, f, cell(size(f, 1), N-1)));
	else
		error('realize() found a field "%s" it cannot handle.', fields{k});
	end
end

if floor(1 / N * 100) > progress
	progress = floor(1 / N * 100);
	fprintf(1, '\b\b\b%02d%%', progress);
end

for r = 2:N
	s = fetch_resource(resource_files{r}, repo);
	for k = 1:length(fields)
		if size(f, 2) == 1
			% eval() seems to be a million times faster than using subsasgn().
			eval(['data.' fields{k} '(:, r) = s.' fields{k} ';']);
		else
			error('realize() found a field "%s" it cannot handle.', fields{k});
		end
	end
	
	if floor(r / N * 100) > progress
		progress = floor(r / N * 100);
		fprintf(1, '\b\b\b%02d%%', progress);
	end
end

fprintf(1, '\n');

data.Meta = qset;
data.Meta = rmfield(data.Meta, 'Resource');
return;






function data = fetch_resource(resource, repo)

data = [];
if resource(1) == '/'
	error 'Bad resource URL.';
end

data = repo.load_cached(resource);
data = data.data;
return;

