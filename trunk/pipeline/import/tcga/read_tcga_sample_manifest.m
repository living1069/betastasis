function manifest = read_tcga_sample_manifest(filepath)

fid = fopen(filepath);
data = textscan(fid, '%*s %*s %*s %*s %s %*s %s', 'HeaderLines', 1);
fclose(fid);

sample_codes = data{1};
files = data{2};
clear data;

manifest = struct('File', {{}}, ...
                  'Sample', {{}});

%tumor_files = {};
%normal_tissue_files = {};
%cellline_files = {};
all_files = {};

seen_samples = containers.Map();

for k = 1:length(sample_codes)
%	if regexp(sample_codes{k}, 'TCGA-\d+-\d+-11')
%		normal_tissue_files{end + 1} = files{k};
%	elseif regexp(sample_codes{k}, 'TCGA-\d+-\d+-20')
%		cellline_files{end + 1} = files{k};
%	elseif regexp(sample_codes{k}, 'TCGA-\d+-\d+-01')
%		if strcmp(files{k}, '0124.CEL'), continue, end
%		tumor_files{end + 1} = files{k};
%	else
%		continue;  % Only add known types of samples to list of all samples.
%	end

	% Skip duplicate samples.
	if seen_samples.isKey(sample_codes{k}), continue, end
	seen_samples(sample_codes{k}) = 1;

	if regexp(sample_codes{k}, 'TCGA-.+')
		manifest.File{end + 1, 1} = files{k};
		manifest.Sample{end + 1, 1} = sample_codes{k};
	end
end

%
%manifest = struct('All', { all_files }, ...
%                  'Tumor', { tumor_files }, ...
%                  'NormalTissue', { normal_tissue_files }, ...
%				  'CellLine', { cellline_files });

