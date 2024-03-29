function meta = meta_tcga_clinical(meta, clin_dir)

clin = read_tcga_clinical_data(clin_dir);
meta = tcga_metadata(meta, clin);









function meta = tcga_metadata(meta, clinical)

if strcmp(meta.sample_id{1}(1:4), 'TCGA')
	patient_ids = meta.sample_id;
	for k = 1:length(patient_ids)
		if length(patient_ids{k}) < 12
			patient_ids{k} = '';
		else
			patient_ids{k} = patient_ids{k}(1:12);
		end
	end
	
	[~, pos] = ismember(patient_ids, clinical.patient_id);
else
	clinical.bcr_patient_uuid(1:10)
	[~, pos] = ismember(meta.sample_id, clinical.bcr_patient_uuid);
end

for k = find(pos == 0)
	fprintf('Records for sample %s were not found.\n', meta.sample_id{k});
end

pos_fake = pos;
pos_fake(pos_fake == 0) = 1;
clinical = filter_rows(clinical, pos_fake);

% Now we empty out all the metadata for missing patients.
fields = fieldnames(clinical);
for k = 1:length(fields)
	if isnumeric(getfield(clinical, fields{k}))
		eval(['clinical.' fields{k} '(pos == 0) = NaN;']);
	elseif iscellstr(getfield(clinical, fields{k}))
		eval(['clinical.' fields{k} '(pos == 0) = ' ...
			'repmat({''''}, 1, sum(pos == 0));']);
	end
end

fields = fieldnames(clinical);
for k = 1:length(fields)
	eval(['meta.' fields{k} ' = clinical.' fields{k} ''';']);
end

for k = 1:length(meta.sample_id)
	meta.sample_type{k} = parse_sample_type(clinical.sample_type{k}, ...
		meta.sample_id{k});
end










function sample_type = parse_sample_type(sample_type, sample_id)

if regexpi(sample_id, 'Stratagene.*Ref')
	sample_type = 'Stratagene universal reference';
	return;
elseif regexpi(sample_id, 'Promega.*Ref')
	sample_type = 'Promega universal reference';
	return;
end

if length(sample_id) < 15 || ~strcmp(sample_id(1:4), 'TCGA')
	sample_type = '-';
	return;
end

typecode = str2double(sample_id(14:15));

typecode_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
typecode_map(1) = 'Solid tumor';
typecode_map(2) = 'Recurrent solid tumor';
typecode_map(3) = 'Primary blood derived cancer';
typecode_map(4) = 'Recurrent blood derived cancer';
typecode_map(6) = 'Metastatic tumor';
typecode_map(10) = 'Normal blood';
typecode_map(11) = 'Adjacent normal tissue';
typecode_map(12) = 'Buccal smear';
typecode_map(13) = 'EBV immortalized normal';
typecode_map(20) = 'Normal cell line';

if ~typecode_map.isKey(typecode)
	fprintf('Unknown tissue type %d in TCGA data.', typecode);
	return;
end

if ~strcmp(sample_type, '-')
	sample_type = [typecode_map(typecode) ', ' sample_type];
else
	sample_type = typecode_map(typecode);
end









function clinical = read_tcga_clinical_data(clinical_dir)

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

fprintf('Reading patient information from file %s...\n', patient_file);
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

clinical = struct;
for k = 1:length(data)
	token = tokens{k}; col_name = token{1};
	if length(col_name) <= 63
		%fprintf(1, 'Adding field %s\n', col_name);
		clinical = setfield(clinical, col_name, replace_nulls(data{k}));
	else
		%fprintf(1, 'Field name too long: %s\n', col_name);
	end
end

S = length(clinical.bcr_patient_barcode);

clinical.patient_id = clinical.bcr_patient_barcode;
clinical = rmfield(clinical, 'bcr_patient_barcode');

capitalize_fields = {'gender', 'vital_status', 'race'};
for k = 1:length(capitalize_fields)
	eval(sprintf('clinical.%s = capitalize(clinical.%s);', ...
		capitalize_fields{k}, capitalize_fields{k}));
end

if isfield(clinical, 'days_to_birth')
	clinical.age_at_diagnosis = -str2double(clinical.days_to_birth) / 365;
	clinical.age_at_death = clinical.age_at_diagnosis + ...
		str2double(clinical.days_to_death) / 365;

	clinical.age_at_last_followup = clinical.age_at_diagnosis + ...
		str2double(clinical.days_to_last_followup) / 365;
	
	if isfield(clinical, 'days_to_tumor_recurrence')
		clinical.age_at_tumor_recurrence = clinical.age_at_diagnosis + ...
			str2double(clinical.days_to_tumor_recurrence) / 365;
	end
	
	if isfield(clinical, 'days_to_tumor_progression')
		clinical.age_at_tumor_progression = clinical.age_at_diagnosis + ...
			str2double(clinical.days_to_tumor_progression) / 365;
	end
end

clinical.survival_time = nan(S, 1);
clinical.survival_time_censored = nan(S, 1);

deceased = strcmpi('deceased', clinical.vital_status);

% Here we could try to rule out patients who died of other causes
% by checking if they died on the day of their last followup, but
% I don't think that's a reliable approach. There's many patients
% for whom death occurs just a couple of days after their last followup.
clinical.survival_time(deceased) = round(365 * ...
	(clinical.age_at_death(deceased) - clinical.age_at_diagnosis(deceased)));
clinical.survival_time_censored(deceased) = 0;

clinical.survival_time(~deceased) = round(365 * ...
	(clinical.age_at_last_followup(~deceased) - ...
	clinical.age_at_diagnosis(~deceased)));
clinical.survival_time_censored(~deceased) = 1;

if isfield(clinical, 'histological_type')
	clinical.sample_type = replace_nulls(clinical.histological_type);
else
	clinical.sample_type = repmat({'-'}, S, 1);
end

if isfield(clinical, 'anatomic_organ_subdivision')
	clinical.organ = replace_nulls(clinical.anatomic_organ_subdivision);
end


	
% Remove obsolete fields that have been converted to a more suitable format.
remove_fields = {'days_to_.*'};
fields = fieldnames(clinical);
should_remove = false(1, length(fields));
for k = 1:length(remove_fields)
	should_remove(rx(fields, remove_fields{k})) = true;
end
for k = find(should_remove)
	clinical = rmfield(clinical, fields{k});
end


patient_id_to_idx = containers.Map(clinical.patient_id, ...
	num2cell(1:length(clinical.patient_id)));


	
	
	
	

	
	


if 0
	
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

clinical.Treatment = repmat({'-'}, length(clinical.ID), 1);

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
	
	if strcmp(clinical.drug_treatment{patient_idx}, '-')
		clinical.drug_treatment{patient_idx} = drug_descr;
	else
		clinical.drug_treatment{patient_idx} = sprintf('%s\n%s', ...
			clinical.drug_treatment{patient_idx}, drug_descr);
	end
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
		cellstr{k} = '';
	end
end

