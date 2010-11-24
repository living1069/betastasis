function [] = check_dataset_path(dataset)

global pipeline_config;
if isempty(pipeline_config.Repositories)
	error 'No repositories have been configured.';
end
root = pipeline_config.Repositories{1}.URL;

if ~exist(root)
	error 'Cannot find pipeline datasets directory.';
end

if exist([root '/' flatten_str(dataset)]) ~= 0
	error 'A dataset with the given name already exists.';
end

% Check for the case where the user is trying to create a dataset within a
% folder, but a dataset already exists with the same name as the folder.
dataset_folder = dataset;
slash_pos = find(dataset == '/');
if ~isempty(slash_pos)
	dataset_folder = dataset(1:slash_pos(end)-1);
end

if exist([root '/' flatten_str(dataset_folder) '/metadata.mat']) ~= 0
	error(['Cannot create a new dataset within that folder, because there ' ...
		'already exists a dataset with the same name as the folder.']);
end

