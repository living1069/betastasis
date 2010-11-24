function metadata = import_tcga_uarray_data(dataset, platform)

samples = {};
array_files = {};
channels = {};

files = dir('.');
for k = 1:length(files)
	fname = files(k).name;
	if fname(1) == '.', continue, end
	
	if files(k).isdir
		[tmp_samples, tmp_array_files, tmp_channels] = ...
			scan_tcga_batch_dir(fname);
		samples = cat(1, samples, tmp_samples);
		array_files = cat(1, array_files, tmp_array_files);
		channels = cat(1, channels, tmp_channels);
	end
end

if length(channels) > 0 && length(channels) < length(samples)
	error 'Some SDRF files are missing channel information.';
end

% Remove duplicate samples and array files.
if ~isempty(channels)
	tmp = strcat(array_files, '(', channels, ')');
	[~, uniq] = unique(tmp);
	samples = samples(uniq);
	array_files = array_files(uniq);
	channels = channels(uniq);
else
	[~, uniq] = unique(array_files);
	samples = samples(uniq);
	array_files = array_files(uniq);
end

fprintf(1, 'Found %d microarray samples in SDRF files.\n', length(array_files));

file_listing = ptemp;
system(['find > ' file_listing]);

all_files = importdata(file_listing);
system(['rm ' file_listing]);

missing = false(length(array_files), 1);

for k = 1:length(array_files)
	offsets = strfind(all_files, array_files{k});
	matches = false(size(offsets));
	for m = 1:length(offsets)
		if ~isempty(offsets{m})
			matches(m) = true;
		end
	end
	
	matches = find(matches);
	if length(matches) > 1
		fprintf(1, 'Array file %s was found in multiple TCGA archives:\n', ...
			array_files{k});
		for m = 1:length(matches)
			fprintf(1, '- %s\n', all_files{matches(m)});
		end
	elseif length(matches) == 0
		fprintf(1, ['Array file %s corresponding to sample %s was not found '...
		            'in any TCGA archive. Skipping it...\n'], ...
			array_files{k}, samples{k});
		missing(k) = true;
		continue;
	end
	
	tmp = all_files{matches(1)};
	tmp = tmp(3:end);   % Remove the ./ prefix
	
	% Remove any possible compressed archive suffix.
	tmp = regexprep(tmp, '\.(bz2|gz)$', '', 'ignorecase');
	
	array_files{k} = tmp;
end

if any(missing)
	samples = samples(~missing);
	array_files = array_files(~missing);
end

sample_ids = samples;
for k = 1:length(sample_ids)
	if length(sample_ids{k}) < 4 || ~strcmp(sample_ids{k}(1:4), 'TCGA')
		continue;
	end
	sample_ids{k} = sample_ids{k}(1:15);
end

if ~isempty(channels)
	file_to_sample = containers.Map( ...
		strcat(array_files, '(', channels, ')'), sample_ids);
else
	file_to_sample = containers.Map(array_files, sample_ids);
end

import_uarray_data(dataset, platform, 'IncludeFiles', unique(array_files), ...
	'SampleMapping', file_to_sample);

metadata = augment_meta_tcga(dataset);

return;





function [samples, array_files, channels] = scan_tcga_batch_dir(dirname)

sdrf = '';
samples = {};
array_files = {};
channels = {};

files = dir(dirname);
for k = 1:length(files)
	if files(k).isdir, continue, end
	
	if regexpi(files(k).name, '.+\.sdrf(\.txt)?$')
		if ~isempty(sdrf)
			error('Found more than one SDRF file in directory %s.', dirname);
		end
		sdrf = files(k).name;
		continue;
	end
end

if ~isempty(sdrf)
	fprintf(1, 'Reading SDRF file %s...\n', [dirname '/' sdrf]);
	[samples, array_files, channels] = read_tcga_sdrf([dirname '/' sdrf]);
end

return;


