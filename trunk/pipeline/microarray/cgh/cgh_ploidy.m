function ploidy = cgh_ploidy(test, ref)

global organism;

if isfield(test, 'meta'), test = test.meta; end
if isfield(ref, 'meta'), ref = ref.meta; end

S = size(test, 2);

if nargin == 2
	genders_match = false(S, 1);
	if isfield(test, 'gender') && isfield(ref, 'gender')
		genders_match = strcmpi(test.gender, ref.gender) & ...
			rx(test.gender, 'Male|Female');
	end

	if ~all(genders_match)
		fprintf(['WARNING: Skipping sex chromosomes in %d samples with ' ...
			'mismatched or missing gender information.\n'], ...
			sum(~genders_match));
	end
elseif nargin == 1
	genders_match = true(S, 1);
end

ploidy = nan(length(organism.Chromosomes.Name), S);
for s = 1:S
	if ~genders_match(s) || strcmp(test.gender{s}, '-')
		ploidy(:, s) = [2*ones(1, 22) NaN NaN NaN]';
	elseif strcmpi(test.gender{s}, 'Male')
		ploidy(:, s) = [2*ones(1, 22) 1 1 NaN]';
	elseif strcmpi(test.gender{s}, 'Female')
		ploidy(:, s) = [2*ones(1, 22) 2 0 NaN]';
	else
		ploidy(:, s) = [2*ones(1, 22) NaN NaN NaN]';
		fprintf(1, 'Ploidy of gender "%s" is unknown.\n', test.gender{s});
	end
end

