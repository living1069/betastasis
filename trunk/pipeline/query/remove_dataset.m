function [] = remove_dataset(name)

global pipeline_config;
repo = pipeline_config.Repositories{1};   % FIXME

matches = find_dataset(name, repo);
if length(matches) == 0
	fprintf(1, 'No dataset by that name exists.\n');
	return;
elseif length(matches) > 1
	fprintf(1, 'Ambiguous dataset name. Possible matches include:\n');
	for k = 1:length(matches)
		fprintf(1, '- %s\n', strrep(matches{k}, '_', ' '));
	end
	return;
end

path = [ppath '/datasets/' matches{1}];

if exist([path '/raw'])
	error(['Dataset contains raw resource files in ' path '/raw. ' ...
	       'Please delete them manually before continuing.']);
end

% We double check here that we really only had exactly one matching dataset.
% "rm -fr" is very dangerous stuff, so a double check is in order in case
% the code above is changed and acquires a bug.
if length(matches) == 1
	system(['rm -fr ' path]);
end

