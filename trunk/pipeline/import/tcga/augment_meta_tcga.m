function metadata = augment_meta_tcga(dataset)

global pipeline_config;

metadata = query(dataset);

tcga_patients = struct;
tcga_samples = struct;
tcga_misc = struct;

files = dir(pipeline_config.TCGA.Path);
for k = 1:length(files)
	if files(k).name(1) == '.' || ~files(k).isdir, continue, end
	
	sources = {'intgen.org', 'nationwidechildrens.org'};
	for s = 1:length(sources)
		clin_dir = [pipeline_config.TCGA.Path '/' files(k).name ... 
			'/bcr/' sources{s} '/biotab/clin'];
		if exist(clin_dir)
			[tmp_patients, tmp_samples, tmp_misc] = ...
				read_tcga_clinical_data(clin_dir);
			tcga_patients = cat_structs(tcga_patients, tmp_patients);
			tcga_samples = cat_structs(tcga_samples, tmp_samples);
			tcga_misc = cat_structs(tcga_misc, tmp_misc);
		end
	end
end

if length(fieldnames(tcga_patients)) == 0
	fprintf(1, 'WARNING: No clinical metadata available.\n');
	return;
end

id_to_index = containers.Map(tcga_patients.ID, ...
	num2cell(1:length(tcga_patients.ID)));
	
patient_ids = metadata.Sample.ID;
for k = 1:length(patient_ids)
	if length(patient_ids{k}) < 12 || ~strcmp(patient_ids{k}(1:4), 'TCGA')
		patient_ids{k} = '-';
	else
		patient_ids{k} = patient_ids{k}(1:12);
	end
end

missing = ~id_to_index.isKey(patient_ids);
indices = zeros(length(patient_ids), 1);

indices(~missing) = cell2mat(id_to_index.values(patient_ids(~missing)));
indices(missing) = NaN;

missing = find(missing);
for k = 1:length(missing)
	if strcmp(patient_ids{missing(k)}, '-'), continue, end
	fprintf(1, 'Records for patient %s were not found.\n', ...
		patient_ids{missing(k)});
end

patients = permute_struct_fields(tcga_patients, indices);
samples = permute_struct_fields(tcga_samples, indices);
misc = permute_struct_fields(tcga_misc, indices);

% The patient IDs for missing patient records need to be returned
% back to their original states.
patients.ID(missing) = patient_ids(missing);
samples.ID = metadata.Sample.ID;

for k = 1:length(samples.ID)
	samples.Type{k} = parse_sample_type(samples.Type{k}, samples.ID{k});
end

samples = associate_subtypes(samples, ...
	[pipeline_config.TCGA.Path '/gbm/tcga_gbm_subtypes.txt']);
samples = associate_subtypes(samples, ...
	[pipeline_config.TCGA.Path '/ov/tcga_ov_subtypes.txt']);

fields = fieldnames(samples);
for k = 1:length(fields)
	eval(['metadata.Sample.' fields{k} ' = samples.' fields{k} ';']);
end

fields = fieldnames(patients);
for k = 1:length(fields)
	if ~isfield(metadata, 'Patient'), metadata.Patient = struct; end
	eval(['metadata.Patient.' fields{k} ' = patients.' fields{k} ';']);
end

fields = fieldnames(misc);
for k = 1:length(fields)
	if ~isfield(metadata, 'Misc'), metadata.Misc = struct; end
	eval(['metadata.Misc.' fields{k} ' = misc.' fields{k} ';']);
end

save_metadata(metadata);








function sample_type = parse_sample_type(sample_type, sample_id)

if regexpi(sample_id, 'Stratagene.*Ref')
	sample_type = 'Stratagene universal reference DNA';
	return;
elseif regexpi(sample_id, 'Promega.*Ref')
	sample_type = 'Promega universal reference DNA';
	return;
end

if length(sample_id) < 15 || ~strcmp(sample_id(1:4), 'TCGA')
	samples_type = '-';
	return;
end

typecode = str2double(sample_id(14:15));
if typecode == 1
	if strcmp(sample_type, '-')
		sample_type = 'Solid tumor';
	else
		sample_type = [sample_type ', solid tumor'];
	end
elseif typecode == 2
	if strcmp(sample_type, '-')
		sample_type = 'Recurrent solid tumor';
	else
		sample_type = [sample_type ', recurrent solid tumor'];
	end
elseif typecode == 3
	if strcmp(sample_type, '-')
		sample_type = 'Primary blood derived cancer';
	else
		sample_type = [sample_type ', primary blood derived cancer'];
	end
elseif typecode == 4
	if strcmp(sample_type, '-')
		sample_type = 'Recurrent blood derived cancer';
	else
		sample_type = [sample_type ', recurrent blood derived cancer'];
	end
elseif typecode == 10
	sample_type = 'Normal blood';
elseif typecode == 11
	if strcmp(sample_type, '-')
		sample_type = 'Adjacent normal tissue';
	else
		sample_type = [sample_type ', adjacent normal tissue'];
	end
elseif typecode == 12
	sample_type = 'Buccal smear';
elseif typecode == 13
	sample_type = 'EBV immortalized normal';
elseif typecode == 20
	sample_type = 'Normal cell line';
else
	fprintf(1, 'Unknown tissue type %d in TCGA data.', typecode);
	sample_type = '-';
end








function samples = associate_subtypes(samples, subtype_file)

fid = fopen(subtype_file);
if fid == -1, return, end
data = textscan(fid, '%s %s');
fclose(fid);

subtype_sample_ids = data{1};
for k = 1:length(subtype_sample_ids)
	if length(subtype_sample_ids{k}) < 12, continue, end
	subtype_sample_ids{k} = subtype_sample_ids{k}(1:12);
end

sample_ids = samples.ID;
for k = 1:length(sample_ids)
	if length(sample_ids{k}) < 12, continue, end
	sample_ids{k} = sample_ids{k}(1:12);
end

sample_to_subtype = containers.Map(subtype_sample_ids, data{2});
found = sample_to_subtype.isKey(sample_ids);

if ~isfield(samples, 'SubtypeTCGA')
	samples.SubtypeTCGA = repmat({'-'}, length(sample_ids), 1);
end
samples.SubtypeTCGA(found) = sample_to_subtype.values(sample_ids(found));


