function [] = update_taylor_patient_metadata()

dirs = dir([ppath '/datasets/taylor_prostate/']);
for k = 1:length(dirs)
	dataset = dirs(k).name;
	if dataset(1) == '.' || ~dirs(k).isdir, continue, end
	
	fprintf(1, 'Updating metadata for dataset %s...\n', ...
		strrep(dataset, '_', ' '));

	qset = query(['taylor_prostate/' dataset]);
	qset.Patient.RecurrenceFreeTime = ...
		30 * str2double(qset.Misc.BiochemicalRecurrence_FreeTime);
	qset.Patient.RecurrenceEvent = ...
		~strcmpi('NO', qset.Misc.BiochemicalRecurrence_Event);
	qset.Patient.Gender = repmat({'Male'}, length(qset.Patient.Gender), 1);
	save_metadata(qset);
end


