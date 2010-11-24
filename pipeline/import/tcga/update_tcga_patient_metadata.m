function [] = update_tcga_patient_metadata()

dirs = dir([ppath '/datasets/tcga/']);
for k = 1:length(dirs)
	dataset = dirs(k).name;
	if dataset(1) == '.' || ~dirs(k).isdir, continue, end
	
	fprintf(1, 'Updating metadata for dataset %s...\n', ...
		strrep(dataset, '_', ' '));

	augment_meta_tcga(['tcga/' dataset]);
end


