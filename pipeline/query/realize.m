function data = realize(qset)

global pipeline_config;

data = struct;

if ~isfield(qset, 'Resource')
	error 'realize() only works on query sets.';
end

N = length(qset.Resource);
if N == 0, error('Query data set is empty.\n'); end

raw_data = false;
if strcmpi('Sequence reads', qset.Type), raw_data = true; end

resources = qset.Resource;
resource_files = cell(size(resources));

% Find a handle to the repository where the resources are stored and
% cache (download) all referred data resources.
repo_ids = {};
for k = 1:length(resources)
	colon = find(resources{k} == ':');
	if length(colon) > 1, error 'Invalid resource identifier.'; end
		
	% Check if we are dealing with a file that is not in any repository.
	if length(colon) == 0
		repo_ids{end+1, 1} = '';
		resource_files{k} = resources{k};
		continue;
	end
		
	repo_ids{end+1, 1} = resources{k}(1:colon-1);
	resource_files{k} = [resources{k}(colon+1:end)];
end

repo_name = unique(repo_ids);
if length(repo_name) > 1
	error 'Realizing multi-repo datasets is not yet supported.';
end
repo_name = repo_name{1};

% FIXME: Kludge. Remove.
if isempty(repo_name) && raw_data
	for k = 1:length(resource_files)
		data.Raw{1, k} = FilePool;
		data.Raw{1, k}.static(resource_files{k});
	end
	data.Meta = rmfield(qset, 'Resource');
	
	return;
end



	
for r = 1:length(pipeline_config.Repositories)
	repo = pipeline_config.Repositories{r};
	if strcmp(repo_name, repo.Name), break, end
	if r == length(pipeline_config.Repositories)
		error('Repository "%s" was not found in configuration.', repo_name);
	end
end




% Raw data resources are handled in a special way.
if raw_data
	raw_groups = [];
	raw_files = {};
	for k = 1:length(resource_files)
		start = find(resource_files{k} == '/');
		if isempty(start)
			start = 1;
			prefix = '';
		else
			start = start(end)+1;
			prefix = resource_files{k}(1:start-1);
		end
		
		ends = [find(resource_files{k} == ';'), length(resource_files{k})+1];
		
		pos = start;
		for f = 1:length(ends)
			raw_groups(end+1, 1) = k;
			raw_files{end+1, 1} = [prefix resource_files{k}(pos:ends(f)-1)];
			pos = ends(f)+1;
		end
	end
	
	cached = repo.cache(raw_files);
	
	progress = Progress;
	
	for k = 1:length(resource_files)
		data.Raw{1, k} = FilePool;
		cache_group = cached(raw_groups == k);
		for f = 1:length(cache_group)
			data.Raw{1, k}.static(cache_group{f});
		end
		progress.update(k / length(resource_files));
	end
	
else
	
	% Normal .mat resources are handled here.
	resource_files{k} = strcat(resource_files{k}, '.mat');
	cached = repo.cache(resource_files);

	progress = Progress;
	
	for d = 1:N
		s = load(cached{d}); s = s.data;
		fields = fieldnames(s);
		
		if d == 1
			% Extend the columns of the first resource into matrices
			% in preparation for loading the other data resources.
			for k = 1:length(fields)
				f = getfield(s, fields{k});
				if isnumeric(f) && size(f, 2) == 1
					eval(['data.' fields{k} ' = f;' ...
						  'data.' fields{k} '(end, N) = 0;']);
				elseif iscell(f) && size(f, 2) == 1
					data = setfield(data, fields{k}, ...
						cat(2, f, cell(size(f, 1), N-1)));
				else
					error('realize() found a field "%s" it cannot handle.', ...
						fields{k});
				end
			end
		else
			for k = 1:length(fields)
				f = getfield(s, fields{k});
				if size(f, 2) == 1
					% eval() is much faster than using subsasgn().
					eval(['data.' fields{k} '(:, d) = s.' fields{k} ';']);
				else
					error('realize() cannot handle field "%s".', fields{k});
				end
			end
		end
	
		progress.update(d/N);
	end
end

data.Meta = rmfield(qset, 'Resource');

