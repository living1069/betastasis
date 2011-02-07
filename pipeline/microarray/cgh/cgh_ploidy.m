function ploidy = cgh_ploidy(test, ref)

global organism;

S = size(test.Mean, 2);

genders_match = false(S, 1);
if isfield(test, 'Meta') && isfield(ref, 'Meta') && ...
	isfield(test.Meta, 'Patient') && isfield(test.Meta.Patient, 'Gender') && ...
	isfield(ref.Meta, 'Patient') && isfield(ref.Meta.Patient, 'Gender')
	
	genders_known = ~strcmpi('-', test.Meta.Patient.Gender) & ...
		~strcmpi('-', ref.Meta.Patient.Gender);
	genders_match = strcmpi(test.Meta.Patient.Gender, ...
		ref.Meta.Patient.Gender) | genders_known;
end

if ~all(genders_match)
	fprintf(1, ['WARNING: Skipping sex chromosomes in %d samples with ' ...
	            'mismatched or missing gender information.\n'], ...
				sum(~genders_match));
end

ploidy = nan(length(organism.Chromosomes.Name), S);
for s = 1:S
	if ~genders_match(s) || strcmp(test.Meta.Patient.Gender{s}, '-')
		ploidy(:, s) = [2*ones(1, 22) NaN NaN NaN]';
	elseif strcmpi(test.Meta.Patient.Gender{s}, 'Male')
		ploidy(:, s) = [2*ones(1, 22) 1 1 NaN]';
	elseif strcmpi(test.Meta.Patient.Gender{s}, 'Female')
		ploidy(:, s) = [2*ones(1, 22) 2 0 NaN]';
	else
		ploidy(:, s) = [2*ones(1, 22) NaN NaN NaN]';
		fprintf(1, 'Ploidy of gender "%s" is unknown.\n', ...
			test.Meta.Patient.Gender{s});
	end
end

