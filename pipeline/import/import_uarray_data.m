function metadata = import_uarray_data(dataset, platform, varargin)

files = {};
sample_map = [];
metafile = '';

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'IncludeFiles')
		files = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'SampleMap')
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
	
	if strcmpi(varargin{k}, 'MetaFile')
		if exist(varargin{k+1}) == 0
			error('Specified clinical metadata file %s does not exist.', ...
				varargin{k+1});
		end
		metafile = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

metadata = struct;
metadata.Type = 'Microarray probe intensities';
metadata.Sample = struct;
metadata.Resource = {};

check_dataset_path(dataset);
	




num_found = 0;

import_method = '';

fprintf(1, 'Loading microarray probe definitions...\n');

if strcmpi(platform, 'Agilent 244K TCGA custom 1')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_244k_tcga_custom_01/probes']);
	probes = probes.agilent_244k_tcga_custom_01_probes;
elseif strcmpi(platform, 'Agilent 244K TCGA custom 2')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_244k_tcga_custom_02/probes']);
	probes = probes.agilent_244k_tcga_custom_02_probes;
elseif strcmpi(platform, 'Agilent 244K TCGA custom 3')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_244k_tcga_custom_03/probes']);
	probes = probes.agilent_244k_tcga_custom_03_probes;
elseif strcmpi(platform, 'Agilent HG CGH 244A')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_hg_cgh_244a/probes']);
	probes = probes.agilent_hg_cgh_244a_probes;
elseif strcmpi(platform, 'Agilent G3 Human CGH 1M')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_g3_cgh_1m/probes']);
	probes = probes.agilent_g3_cgh_1m_probes;
elseif strcmpi(platform, 'Agilent Human miRNA 8x15K')
	import_method = 'agilent_single';
	probes = load([ppath '/platforms/agilent_human_mirna_8x15k/probes']);
	probes = probes.agilent_human_mirna_8x15k_probes;
elseif strcmpi(platform, 'Agilent Human miRNA 8x15K v2')
	import_method = 'agilent_single';
	probes = load([ppath '/platforms/agilent_human_mirna_8x15k_v2/probes']);
	probes = probes.agilent_human_mirna_8x15k_v2_probes;
elseif strcmpi(platform, 'Affymetrix HG U133A')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_hg_u133a/probes']);
	probes = probes.affy_hg_u133a_probes;
elseif strcmpi(platform, 'Affymetrix U95A v2')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_u95a_v2/probes']);
	probes = probes.probes;
elseif strcmpi(platform, 'Affymetrix U95B')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_u95b/probes']);
	probes = probes.probes;
elseif strcmpi(platform, 'Affymetrix U95C')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_u95c/probes']);
	probes = probes.probes;
elseif strcmpi(platform, 'Affymetrix HT HG U133A')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_ht_hg_u133a/probes']);
	probes = probes.affy_ht_hg_u133a_probes;
elseif strcmpi(platform, 'Affymetrix HG U133 Plus 2.0')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_hg_u133_plus_2/probes']);
	probes = probes.affy_hg_u133_plus_2_probes;
elseif strcmpi(platform, 'Affymetrix Human Exon 1.0 ST')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_huex_1_0_st/probes']);
	probes = probes.affy_huex_1_0_st_probes;
elseif strcmpi(platform, 'Affymetrix GW SNP 6')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_gw_snp_6/probes']);
	probes = probes.affy_gw_snp_6_probes;
elseif strcmpi(platform, 'Affymetrix Human Mapping 250K Sty')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_mapping_250k_sty/probes']);
	probes = probes.probes;
elseif strcmpi(platform, 'Affymetrix Human Mapping 250K Nsp')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_mapping_250k_nsp/probes']);
	probes = probes.probes;
else
	error 'Unrecognized array platform specified.';
end

uarray_regexp = '^(.+)\.(cel|txt)(\.bz2|\.gz)?$';
compressed_regexp = '^(.+)(\.bz2|\.gz)?$';

if isempty(files)
	entries = dir('.');
	for k = 1:length(entries)
		if entries(k).isdir, continue, end
	
		tokens = regexpi(entries(k).name, uarray_regexp, 'tokens');
		if isempty(tokens), continue, end
			
		files{end + 1, 1} = entries(k).name;
	end
else
	% If we have a list of user specified filenames, we need to check if the
	% files exist in a compressed form.
	for k = 1:length(files)
		if exist(files{k}), continue, end
		
		if exist([files{k}, '.bz2'])
			files{k} = [files{k} '.bz2'];
		elseif exist([files{k}, '.gz'])
			files{k} = [files{k} '.gz'];
		end
	end
	
	files = files(:);
end

if length(files) == 0
	fprintf(1, 'No microarray samples found.\n');
	return;
end

N_samples = 0;
metadata.Sample.Filename = {};
if strcmp(import_method, 'agilent_dual')
	metadata.Sample.Channel = {};
end

for k = 1:length(files)
	% Extract the compressed sample file if necessary.
	file = files{k};
	compressed = false;
	
	if exist(files{k}) == 0
		fprintf(1, 'File %s does not exist. Skipping it...\n', files{k});
		failed(k) = true;
		continue;
	end
	
	tokens = regexpi(files{k}, '^(.+)\.gz$', 'tokens');
	if length(tokens) == 1
		tokens = tokens{1}; file = tokens{1};
		compressed = true;
		system(['gunzip -c ' files{k} ' > ' file]);
	end
	
	tokens = regexpi(files{k}, '^(.+)\.bz2$', 'tokens');
	if length(tokens) == 1
		tokens = tokens{1}; file = tokens{1};
		compressed = true;
		system(['bunzip2 -k ' files{k}]);
	end

	fprintf(1, 'Importing %s...\n', file);
	
	try
		data = struct;
		if strcmp(import_method, 'agilent_dual')
			channels = read_uarray_sample_agilent(file, probes, ...
				{'gProcessedSignal', 'rProcessedSignal'});
			data.GreenMean = channels(:, 1);
			data.RedMean = channels(:, 2);
		elseif strcmp(import_method, 'agilent_single')
			data.Mean = read_uarray_sample_agilent(file, probes);
		elseif strcmp(import_method, 'affy_cel')
			data = read_uarray_sample_cel(file, probes);
		end
	catch exception
		fprintf(1, 'Skipping sample due to import failure. Reason:\n');
		fprintf(1, '%s\n', exception.message);
		continue;
	end
	
	if compressed
		system(['rm ' file]);
	end
	
	if strcmp(import_method, 'agilent_dual')
		dual = data;
		
		N_samples = N_samples + 1;
		data = struct; data.Mean = dual.GreenMean;
		metadata.Resource{N_samples, 1} = create_sample(dataset, data);
		metadata.Sample.Filename{N_samples, 1} = file;
		metadata.Sample.Channel{N_samples, 1} = 'Cy3';
		
		N_samples = N_samples + 1;
		data = struct; data.Mean = dual.RedMean;
		metadata.Resource{N_samples, 1} = create_sample(dataset, data);
		metadata.Sample.Filename{N_samples, 1} = file;
		metadata.Sample.Channel{N_samples, 1} = 'Cy5';
	else
		N_samples = N_samples + 1;
		metadata.Resource{N_samples, 1} = create_sample(dataset, data);
		metadata.Sample.Filename{N_samples, 1} = file;
	end
end

if ~isempty(sample_map)
	if strcmp(import_method, 'agilent_dual')
		keys = strcat(metadata.Sample.Filename, ...
			'(', metadata.Sample.Channel, ')');
	else
		keys = metadata.Sample.Filename;
	end
	
	found = sample_map.isKey(keys);
	missing = find(~found);
	for k = 1:length(missing)
		fprintf(1, 'No sample ID association for "%s".\n', keys{missing(k)});
	end
	
	metadata.Sample.ID = repmat({'-'}, length(keys), 1);
	metadata.Sample.ID(found) = sample_map.values(keys(found));
end

metadata.Platform = repmat({platform}, length(metadata.Resource), 1);

save([ppath '/datasets/' flatten_str(dataset) '/metadata.mat'], 'metadata');

if ~isempty(metafile)
	metadata = augment_meta_tab(dataset, metafile);
end

