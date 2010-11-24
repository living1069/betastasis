
% QUERY       Query dataset repositories for samples fulfilling some criteria.
%
%    QSET = QUERY(DS) searches all configured dataset repositories for a
%    dataset with the name DS. The dataset name can be specified in full, or
%    only partially. Use of wildcards is allowed in the dataset name DS. A
%    query set referencing all samples from the dataset is retrieved. 
%    If the search string DS does not uniquely identify a dataset, then
%    [] is returned instead and a warning is shown.
%
%    QSET = QUERY(DS, PRED) retrieves a dataset, and then filters the samples
%    according to the predicate PRED. The predicate can be specified as a
%    query language string (i.e. 'sample type = primary tumor'), logical
%    vector or index vector. A query set referencing the samples that fulfill
%    the predicate is returned.
%   
%    See also FILTER_QUERY for a list of available predicate keys.

% Author: Matti Annala <matti.annala@tut.fi>

function qset = query(dataset, predicate)

global pipeline_config;
if isempty(pipeline_config.Repositories)
	error 'No repositories have been configured.';
end

if nargin == 0, dataset = ''; end
if nargin < 2, predicate = ''; end
	
total_repos = {};
total_items = {};

for r = 1:length(pipeline_config.Repositories)
	repo = pipeline_config.Repositories{r};
	items = find_dataset(dataset, repo);

	if isempty(dataset) && ~isempty(items)
		fprintf(1, 'Available data sets in %s repository:\n', repo.Name);

		for k = 1:length(items)
			slashes = find(items{k} == '/');
			if isempty(slashes) || ...
				(length(slashes) == 1 && slashes(1) == length(items{k}))
				fprintf(1, '- %s\n', strrep(items{k}, '_', ' '));
			end
		end
		
		fprintf(1, '\n');
		continue;
	end
	
	total_repos = cat(1, total_repos, repmat({repo.Name}, length(items), 1));
	total_items = cat(1, total_items, items);
end

if isempty(dataset), return, end

if length(total_items) == 0
	error 'No datasets by that name could be found.';
elseif length(total_items) > 1
	fprintf(1, 'Ambiguous query.\n');
	for r = 1:length(pipeline_config.Repositories)
		repo = pipeline_config.Repositories{r};
		repo_items = total_items(strcmp(repo.Name, total_repos));
		if isempty(repo_items), continue, end
		
		fprintf(1, 'Possible matches in %s repository:\n', repo.Name);
		for k = 1:length(repo_items)
			fprintf(1, '- %s\n', strrep(repo_items{k}, '_', ' '));
		end
		fprintf(1, '\n');
	end

	return;
end

% If we arrive here, then we have found a unique match for the dataset query.
% Now we find a handle to the repository where the matching dataset was found.
for r = 1:length(pipeline_config.Repositories)
	repo = pipeline_config.Repositories{r};
	if strcmp(repo.Name, total_repos{1})
		ds = total_items{1};
		metadata = repo.load([ds '/metadata.mat']);
		qset = metadata.metadata;
		
		% Check if the dataset is using old style resource URLs.
		res = qset.Resource{1};
		if any(res == '/') || any(res == ':')
			%fprintf(1, ['WARNING: Dataset uses old style resource URLs. ' ...
			%	'Converting to new style URLs...\n']);
			qset.Resource = regexprep(qset.Resource, '^.+/(.+?)$', '$1');
		end

		
		% Append the repository name to the resource identifiers.
		qset.Resource = strcat(repo.Name, ':', ds, '/', qset.Resource);
		break;
	end
end

qset = filter_query(qset, predicate);

