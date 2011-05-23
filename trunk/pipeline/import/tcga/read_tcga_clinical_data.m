function [patients, samples, misc] = read_tcga_clinical_data(clinical_dir)

files = dir(clinical_dir);

patient_file = '';
drug_file = '';

for k = 1:length(files)
	if regexp(files(k).name, 'patient_all')
		patient_file = files(k).name;
	end
	
	if regexp(files(k).name, 'clinical_drug_all')
		drug_file = files(k).name;
	end
end

if isempty(patient_file)
	error('No patient file found in %s', clinical_dir);
end

fprintf(1, 'Reading patient information from file %s...\n', patient_file);
fid = fopen([clinical_dir '/' patient_file]);

line = fgetl(fid);
if line == -1, return, end
	
% Figure out which columns we should read.
tokens = regexp(line, '(\w+)(\t|$)', 'tokens');

% Construct a format string for parsing.
parse_format = '';
for k = 1:length(tokens)
	parse_format = [parse_format '%s'];
end

data = textscan(fid, parse_format, 'Delimiter', '\t');
fclose(fid);

misc = struct;
for k = 1:length(data)
	token = tokens{k}; col_name = token{1};
	if length(col_name) <= 63
		%fprintf(1, 'Adding field %s\n', col_name);
		misc = setfield(misc, col_name, replace_nulls(data{k}));
	else
		%fprintf(1, 'Field name too long: %s\n', col_name);
	end
end

S = length(misc.bcr_patient_barcode);

patients = struct;
patients.ID = misc.bcr_patient_barcode;

if isfield(misc, 'gender')
	patients.Gender = capitalize(misc.gender);
end

if isfield(misc, 'vital_status')
	patients.Status = capitalize(misc.vital_status);
end

if isfield(misc, 'race')
	patients.Race = capitalize(misc.race);
end

patients.AgeAtDiagnosis = nan(length(patients.ID), 1);
patients.AgeAtDeath = nan(length(patients.ID), 1);
patients.AgeAtLastFollowup = nan(length(patients.ID), 1);
patients.AgeAtTumorRecurrence = nan(length(patients.ID), 1);
patients.AgeAtTumorProgression = nan(length(patients.ID), 1);

for k = 1:length(patients.ID)
	if isfield(misc, 'days_to_birth')
		patients.AgeAtDiagnosis(k) = -str2double(misc.days_to_birth{k}) / 365;
	end
	
	if isfield(misc, 'days_to_death')
		patients.AgeAtDeath(k) = patients.AgeAtDiagnosis(k) + ...
			str2double(misc.days_to_death{k}) / 365;
	end
	
	if isfield(misc, 'days_to_last_followup')
		patients.AgeAtLastFollowup(k) = patients.AgeAtDiagnosis(k) + ...
			str2double(misc.days_to_last_followup{k}) / 365;
	end
	
	if isfield(misc, 'days_to_tumor_recurrence')
		patients.AgeAtTumorRecurrence(k) = patients.AgeAtDiagnosis(k) + ...
			str2double(misc.days_to_tumor_recurrence{k}) / 365;
	end
	
	if isfield(misc, 'days_to_tumor_progression')
		patients.AgeAtTumorProgression(k) = patients.AgeAtDiagnosis(k) + ...
			str2double(misc.days_to_tumor_progression{k}) / 365;
	end
end

if isfield(patients, 'Status')
	patients.SurvivalTime = nan(length(patients.ID), 1);
	patients.Censored = nan(length(patients.ID), 1);
	for k = 1:length(patients.SurvivalTime)
		if strcmp('Deceased', patients.Status{k})
			% Here we could try to rule out patients who died of other causes
			% by checking if they died on the day of their last followup, but
			% I don't think that's a reliable approach. There's many patients
			% for whom death occurs just a couple of days after their last
			% followup.
			patients.SurvivalTime(k) = round(365 * ...
				(patients.AgeAtDeath(k) - patients.AgeAtDiagnosis(k)));
			patients.Censored(k) = 0;
		else
			patients.SurvivalTime(k) = round(365 * ...
				(patients.AgeAtLastFollowup(k) - patients.AgeAtDiagnosis(k)));
			patients.Censored(k) = 1;
		end
	end
end

samples = struct;

if isfield(misc, 'histological_type')
	samples.Type = replace_nulls(misc.histological_type);
else
	samples.Type = repmat({'-'}, S, 1);
end

if isfield(misc, 'anatomic_organ_subdivision')
	samples.Organ = replace_nulls(misc.anatomic_organ_subdivision);
end

samples.Type = regexprep(samples.Type, '\(De Nova\) ', '');

if isfield(patients, 'Race')
	patients.Race = regexprep(patients.Race, ' or african american', '');
end

patient_id_to_idx = containers.Map(patients.ID, ...
	num2cell(1:length(patients.ID)));




	
	
if isempty(drug_file)
	fprintf(1, 'No drug file found in %s...\n', clinical_dir);
	return;
end
	
fprintf(1, 'Reading drug information from file %s...\n', drug_file);
fid = fopen([clinical_dir '/' drug_file]);

line = fgetl(fid);
if line == -1, return, end
	
% Figure out which columns we should read.
tokens = regexp(line, '(\w+)(\t|$)', 'tokens');

% Construct a format string for parsing.
parse_format = '';
for k = 1:length(tokens)
	parse_format = [parse_format '%s'];
end

data = textscan(fid, parse_format, 'Delimiter', '\t');
fclose(fid);

drugs = struct;
for k = 1:length(data)
	token = tokens{k}; col_name = token{1};
	drugs = setfield(drugs, col_name, replace_nulls(data{k}));
end

for k = 1:length(drugs.drug_name)
	str = drugs.drug_name{k};
	if lower(str(1)) == str(1)
		str(1) = upper(str(1));
		drugs.drug_name{k} = str;
	end
end

patients.Treatment = repmat({'-'}, length(patients.ID), 1);

for k = 1:length(drugs.bcr_patient_barcode)
	if ~patient_id_to_idx.isKey(drugs.bcr_patient_barcode{k}), continue, end
		
	patient_idx = patient_id_to_idx(drugs.bcr_patient_barcode{k});
	
	drug_cat = drugs.drug_category{k};
	if strcmp(drug_cat, 'CHEMO'), drug_cat = 'Chemotherapy';
	elseif strcmp(drug_cat, 'HORMONAL'), drug_cat = 'Hormonal treatment';
	elseif strcmp(drug_cat, 'IMMUNO'), drug_cat = 'Immunological treatment';
	elseif strcmp(drug_cat, 'TARGETED'), drug_cat = 'Targeted drug treatment';
	end
	
	if strcmp(drug_cat, '-'), continue, end
		
	if strcmp(drugs.drug_name{k}, '-')
		continue;
	end
		
	drug_time = '';
	%if ~strcmp(drugs.DAYOFDRUGTREATMENTSTART{k}, '-') && ...
	%	~strcmp(drugs.DAYOFDRUGTREATMENTEND{k}, '-')
	%	drug_time = sprintf(', %s - %s', datestr(datenum( ...
	%		str2double(drugs.YEAROFDRUGTREATMENTSTART{k}), ...
	%		str2double(drugs.MONTHOFDRUGTREATMENTSTART{k}), ...
	%		str2double(drugs.DAYOFDRUGTREATMENTSTART{k})), 'yyyy-mm-dd'), ...
	%		datestr(datenum( ...
	%		str2double(drugs.YEAROFDRUGTREATMENTEND{k}), ...
	%		str2double(drugs.MONTHOFDRUGTREATMENTEND{k}), ...
	%		str2double(drugs.DAYOFDRUGTREATMENTEND{k})), 'yyyy-mm-dd'));
	%end
	
	drug_dosage = '';
	if ~strcmp(drugs.drug_dosage{k}, '-')
		if ~strcmp(drugs.dosage_units{k}, '-')
			drug_dosage = sprintf(' %s %s', drugs.drug_dosage{k}, ...
				drugs.dosage_units{k});
		%else
		%	drug_dosage = sprintf(' %s', drugs.drug_dosage{k});
		end
	end

	drug_descr = sprintf('- %s: %s%s%s', drug_cat, drugs.drug_name{k}, ... 
		drug_dosage, drug_time);
	
	if strcmp(patients.Treatment{patient_idx}, '-')
		patients.Treatment{patient_idx} = drug_descr;
	else
		patients.Treatment{patient_idx} = sprintf('%s\n%s', ...
			patients.Treatment{patient_idx}, drug_descr);
	end
end






function cellstr = capitalize(cellstr)
for k = 1:length(cellstr)
	str = lower(cellstr{k});
	str(1) = upper(str(1));
	cellstr{k} = str;
end



function cellstr = replace_nulls(cellstr)
for k = 1:length(cellstr)
	if strcmp(cellstr{k}, 'null')
		cellstr{k} = '-';
	end
end

