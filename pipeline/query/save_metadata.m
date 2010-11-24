function [] = save_metadata(metadata)

if ~isfield(metadata, 'Resource')
	error 'The first argument to save_metadata() must be a query set.';
end

sources = cell(length(metadata.Resource), 1);
for k = 1:length(metadata.Resource)
	tokens = regexp(metadata.Resource{k}, '^(.+)/(.+?)$', 'tokens');
	tokens = tokens{1};
	sources{k} = tokens{1};
	metadata.Resource{k} = tokens{2};
end

if length(unique(sources)) > 1
	error(['save_metadata() does not work on query sets merged from ' ...
		'multiple datasets.']);
end

source = sources{1};
tokens = regexp(source, '^(.+):(.+?)$', 'tokens'); tokens = tokens{1};
dataset = tokens{2};

path = [ppath '/datasets/' dataset];
if ~exist(path), mkdir(path); end

save([path '/metadata.mat'], 'metadata');

