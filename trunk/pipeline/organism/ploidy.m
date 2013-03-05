function ploidy = ploidy(gender, ref_gender)

global organism;

S = length(gender);

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
	if ~genders_match(s) || strcmp(gender{s}, '-')
		ploidy(:, s) = [2*ones(1, 22) NaN NaN NaN]';
	elseif rx(gender{s}, '^male')
		ploidy(:, s) = [2*ones(1, 22) 1 1 NaN]';
	elseif rx(gender{s}, 'female')
		ploidy(:, s) = [2*ones(1, 22) 2 0 NaN]';
	else
		ploidy(:, s) = [2*ones(1, 22) NaN NaN NaN]';
		fprintf('Ploidy of gender "%s" is unknown.\n', gender{s});
	end
end

