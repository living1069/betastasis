function matches = find_dataset(dataset, repo)

items = sort(repo.contents);

% User did not specify a dataset name, so we return them all.
if nargin == 0 || isempty(dataset)
	matches = items;
	return;
end

dataset = flatten_str(dataset);
dataset = regexptranslate('wildcard', dataset);
if nargin < 2, predicate = ''; end

matches = {};
for k = 1:length(items)
	if regexpi(items{k}, dataset)
		matches{end + 1, 1} = items{k};
	end
end

if length(matches) > 1
	% If we have multiple matches, try to resolve the ambiguity.
	exact_match = find(strcmpi(dataset, matches));
	if length(exact_match) == 1 && matches{exact_match(1)}(end) ~= '/'
		matches = matches(exact_match(1));
	end
end

