%
% AUGMENT_META_TAB    Augment dataset with metadata from a tab delimited file.
%    
%    META = AUGMENT_META_TAB(DS, TLM) reads patient and sample specific metadata
%    from the tab delimited file TLM, and associates the metadata with
%    corresponding samples in the pipeline dataset DS.
%    
%    The tab delimited file TLM must begin with a header line that specifies
%    the contents of each column in the file. Supported column names include:
%    
%      General
%      =======
%      FILENAME:           Name of the original sample file.
%      SAMPLE ID:          Sample identifier, shared by technical replicates.
%      REPLICATE ID:       Replicate identifier, unique within a sample.
%      SAMPLE TYPE:        Sample type (tumor/normal/reference DNA etc).
%      SAMPLE SOURCE:      Organization or individual that extracted the sample.
%      TISSUE:             Tissue or organ from which the sample originated.
%      PATIENT ID:         Patient identifier.
%      GENDER:             Patient gender (at birth).
%      RACE:               Patient race.
%      CLINICAL STAGE:     Clinical cancer stage.
%      PATHOLOGICAL STAGE: Pathological cancer stage.
%      AGE AT DEATH:       Patient age at death (in years).
%      AGE AT DX:          Patient age at initial diagnosis (in years).
%      AGE AT LF:          Patient age at the last followup (in years).
%      SURVIVAL TIME:      Time between initial diagnosis and time of death or
%                          last followup, depending on whether the patient died
%                          during the experiment.
%      CENSORED:           Whether the survival time is right-censored or not.
%
%      Cancer specific
%      ===============
%      PRIMARY GG:         Primary Gleason grade.
%      SECONDARY GG:       Secondary Gleason grade.
%      TERTIARY GG:        Tertiary Gleason grade.
%      TOTAL GG:           Total Gleason grade (primary + secondary).
%      BIOPSY GG1:         Primary Gleason grade based on biopsy.
%      BIOPSY GG2:         Secondary Gleason grade based on biopsy.
%      BIOPSY GG3:         Tertiary Gleason grade based on biopsy.
%      BIOPSY GGS:         Total Gleason grade based on biopsy.
%      PATHOLOGICAL GG1:   Primary Gleason grade based on pathol. evidence.
%      PATHOLOGICAL GG2:   Secondary Gleason grade based on pathol. evidence.
%      PATHOLOGICAL GG3:   Tertiary Gleason grade based on pathol. evidence.
%      PATHOLOGICAL GGS:   Total Gleason grade based on pathol. evidence.
% 

function metadata = augment_meta_tab(dataset, filepath, varargin)

sample_map = [];

for k = 1:2:length(varargin)
	if strcmpi('SampleMap', varargin{k})
		if isa(varargin{k+1}, 'containers.Map')
			sample_map = varargin{k+1};
		elseif ischar(varargin{k+1})
			error 'Loading of sample mappings from files not supported yet.';
		else
			error(['Sample mapping must be provided as either a ' ...
				'containers.Map or a filename to a tab delimited file.']);
		end
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

metadata = query(dataset);

if ~isempty(sample_map)
	keys = metadata.Sample.Filename;
	if isfield(metadata.Sample, 'Channel')
		keys = strcat(keys, '(', metadata.Sample.Channel, ')');
	end

	found = sample_map.isKey(keys);
	missing = find(~found);
	for k = 1:length(missing)
		fprintf(1, 'No sample ID association for "%s".\n', keys{missing(k)});
	end
	
	metadata.Sample.ID = repmat({'-'}, length(keys), 1);
	metadata.Sample.ID(found) = sample_map.values(keys(found));
end

fid = fopen(filepath);

header = fgetl(fid);
if header == -1, error 'Empty file.'; end

column_names = textscan(header, '%s', -1, 'Delimiter', '\t');
column_names = column_names{1};

parse_format = '';
for k = 1:length(column_names)
	parse_format = [parse_format '%s'];
end

data = textscan(fid, parse_format, 'Delimiter', '\t');
fclose(fid);



% Check if the metadata file contains a resource ID or sample ID column.
filename_regexp = '^File.?name$|Array.?File|Sequence.?File';
sample_id_regexp = '^Sample.?ID$|^Sample.?Name$';

filename_col = 0;
sample_id_col = 0;
	
for k = 1:length(column_names)
	if regexpi(column_names{k}, filename_regexp)
		if filename_col, error 'Two filename columns found.'; end
		filename_col = k;
	end
	
	if regexpi(column_names{k}, sample_id_regexp)
		if sample_id_col, error 'Two sample ID columns found.'; end
		sample_id_col = k;
	end
end

if filename_col == 0 && sample_id_col == 0
	error 'Neither a filename nor sample ID column was found.';
end



% Construct a mapping between the dataset resources and metadata entries.
if filename_col && isfield(metadata.Sample, 'Filename')
	aug_filenames = metadata.Sample.Filename;
	
	if isfield(metadata.Sample, 'Channel')
		aug_filenames = strcat(aug_filenames, '(', metadata.Sample.Channel, ...
			')');
	end
	
	N = length(aug_filenames);
	sheet_key_to_idx = containers.Map(data{filename_col}, ...
		num2cell(1:length(data{filename_col})));
	found = sheet_key_to_idx.isKey(aug_filenames);
	idx = cell2mat(sheet_key_to_idx.values(aug_filenames(found)));
elseif sample_id_col && isfield(metadata.Sample, 'ID')
	N = length(metadata.Sample.ID);
	sheet_key_to_idx = containers.Map(data{sample_id_col}, ...
		num2cell(1:length(data{sample_id_col})));
	found = sheet_key_to_idx.isKey(metadata.Sample.ID);
	idx = cell2mat(sheet_key_to_idx.values(metadata.Sample.ID(found)));
end



% Remove metadata entries that do not correspond to any dataset items.
for k = 1:length(column_names)
	data{k} = data{k}(idx);
end



if ~isfield(metadata, 'Patient'), metadata.Patient = struct; end
if ~isfield(metadata, 'Sample'), metadata.Sample = struct; end

for k = 1:length(column_names)
	if regexpi(column_names{k}, filename_regexp), continue, end
	
	if regexpi(column_names{k}, sample_id_regexp)
		metadata.Sample.ID = repmat({'-'}, N, 1);
		metadata.Sample.ID(found) = data{k};
	end
	
	if regexpi(column_names{k}, 'Patient.?ID')
		metadata.Patient.ID = repmat({'-'}, N, 1);
		metadata.Patient.ID(found) = data{k};
		continue;
	end
	
	if regexpi(column_names{k}, '^(Patient.?)?Age.?At.?Death$')
		metadata.Patient.AgeAtDeath = nan(N, 1);
		metadata.Patient.AgeAtDeath(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, '^(Patient.?)?Age.?At.?(Diagnosis|DX)$')
		metadata.Patient.AgeAtDiagnosis = nan(N, 1);
		metadata.Patient.AgeAtDiagnosis(found) = str2double(data{k});
		continue;
	end

	if regexpi(column_names{k}, '^(Patient.?)?Age.?At.?(Last.?Followup|LF)$')
		metadata.Patient.AgeAtLastFollowup = nan(N, 1);
		metadata.Patient.AgeAtLastFollowup(found) = str2double(data{k});
		continue;
	end

	if regexpi(column_names{k}, '(Patient.?)?Survival.?Time')
		metadata.Patient.SurvivalTime = nan(N, 1);
		metadata.Patient.SurvivalTime(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, '(Patient.?)?Gender')
		gender = data{k};
		for k = 1:length(gender)
			if regexpi(gender{k}, 'Male|M|Man|XY')
				gender{k} = 'Male';
			elseif regexpi(gender{k}, 'Female|F|Woman|XX')
				gender{k} = 'Female';
			else
				gender{k} = '-';
			end
		end
		metadata.Patient.Gender = repmat({'-'}, N, 1);
		metadata.Patient.Gender(found) = gender;
		continue;
	end
	
	if regexpi(column_names{k}, 'Sample.?Type')
		metadata.Sample.Type = repmat({'-'}, N, 1);
		metadata.Sample.Type(found) = beautify(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, '(Sample.?)?Metastasis.?Site')
		metadata.Sample.MetastasisSite = repmat({'-'}, N, 1);
		metadata.Sample.MetastasisSite(found) = beautify(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, '(Patient.?)?Race')
		metadata.Patient.Race = repmat({'-'}, N, 1);
		metadata.Patient.Race(found) = beautify(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, '(Patient.?)?Censored')
		metadata.Patient.Censored = nan(N, 1);
		metadata.Patient.Censored(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, 'Sample.?Organ|Tissue|Tissue.?Type|Organ')
		metadata.Sample.Organ = repmat({'-'}, N, 1);
		metadata.Sample.Organ(found) = beautify(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, '(Sample.?)?Source')
		metadata.Sample.Source = repmat({'-'}, N, 1);
		metadata.Sample.Source(found) = beautify(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, 'Clinical.?Stage')
		metadata.Sample.ClinicalStage = repmat({'-'}, N, 1);
		metadata.Sample.ClinicalStage(found) = beautify(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, 'Pathological.?Stage')
		metadata.Sample.PathologicalStage = repmat({'-'}, N, 1);
		metadata.Sample.PathologicalStage(found) = beautify(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, ...
		'BiopsyPrimaryGleason(Grade)?|BiopsyPrimaryGG|BiopsyGG1|BiopsyG1')
		metadata.Sample.BiopsyPrimaryGleason = nan(N, 1);
		metadata.Sample.BiopsyPrimaryGleason(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, ...
		'BiopsySecondaryGleason(Grade)?|BiopsySecondaryGG|BiopsyGG2|BiopsyG2')
		metadata.Sample.BiopsySecondaryGleason = nan(N, 1);
		metadata.Sample.BiopsySecondaryGleason(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, ...
		'BiopsyTertiaryGleason(Grade)?|BiopsyTertiaryGG|BiopsyGG3|BiopsyG3')
		metadata.Sample.BiopsyTertiaryGleason = nan(N, 1);
		metadata.Sample.BiopsyTertiaryGleason(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, ...
		'BiopsyTotalGleason(Grade)?|BiopsyTotalGG|BiopsyGGT|BiopsyGGS|BiopsyGS')
		metadata.Sample.BiopsyTotalGleason = nan(N, 1);
		metadata.Sample.BiopsyTotalGleason(found) = str2double(data{k});
		continue;
	end

	if regexpi(column_names{k}, ['Path(ological)?PrimaryGleason(Grade)?|' ...
		'Path(ological)?PrimaryGG|Path(ological)?GG1|Path(ological)?G1'])
		metadata.Sample.PathologicalPrimaryGleason = nan(N, 1);
		metadata.Sample.PathologicalPrimaryGleason(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, ['Path(ological)?SecondaryGleason(Grade)?|' ...
		'Path(ological)?SecondaryGG|Path(ological)?GG2|Path(ological)?G2'])
		metadata.Sample.PathologicalSecondaryGleason = nan(N, 1);
		metadata.Sample.PathologicalSecondaryGleason(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, ['Path(ological)?TertiaryGleason(Grade)?|' ...
		'Path(ological)?TertiaryGG|Path(ological)?GG3|Path(ological)?G3'])
		metadata.Sample.PathologicalTertiaryGleason = nan(N, 1);
		metadata.Sample.PathologicalTertiaryGleason(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, ['Path(ological)?TotalGleason(Grade)?|' ...
		'Path(ological)?TotalGG|Path(ological)?GGT|Path(ological)?GGS|' ...
		'Path(ological)?GS'])
		metadata.Sample.PathologicalTotalGleason = nan(N, 1);
		metadata.Sample.PathologicalTotalGleason(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, 'PrimaryGleason(Grade)?|PrimaryGG|GG1|G1')
		metadata.Sample.PrimaryGleason = nan(N, 1);
		metadata.Sample.PrimaryGleason(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, 'SecondaryGleason(Grade)?|SecondaryGG|GG2|G2')
		metadata.Sample.SecondaryGleason = nan(N, 1);
		metadata.Sample.SecondaryGleason(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, 'TertiaryGleason(Grade)?|TertiaryGG|GG3|G3')
		metadata.Sample.TertiaryGleason = nan(N, 1);
		metadata.Sample.TertiaryGleason(found) = str2double(data{k});
		continue;
	end
	
	if regexpi(column_names{k}, 'TotalGleason(Grade)?|TotalGG|GGT|GGS|GS')
		metadata.Sample.TotalGleason = nan(N, 1);
		metadata.Sample.TotalGleason(found) = str2double(data{k});
		continue;
	end

	field = escape_fieldname(column_names{k});
	if strcmp(field, ''), continue, end
		
	if ~isfield(metadata, 'Misc'), metadata.Misc = struct; end
	eval(['metadata.Misc.' field ' = repmat({''-''}, N, 1);']);
	eval(['metadata.Misc.' field '(found) = beautify(data{k});']);
end

save_metadata(metadata);
return;



function ret = beautify(c)
for k = 1:length(c)
	if isempty(c{k})
		c{k} = '-';
		continue;
	end
	if regexpi(c{k}, '^(UNKNOWN|NA|N/A)$')
		c{k} = '-';
		continue;
	end
end
ret = regexprep(c, '^"(.*)"$', '$1');
return;



function ret = escape_fieldname(str)
ret = strrep(str, ' ', '_');
ret = strrep(str, '-', '_');
return;

