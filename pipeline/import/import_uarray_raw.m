
% IMPORT_UARRAY_DATA   Import raw microarray data into the pipeline
%
%    METADATA = IMPORT_UARRAY_RAW(PLATFORM) searches the current
%    working directory for raw microarray sample files, and imports any found
%    samples into the pipeline as a dataset. The name of the new dataset is
%    specified by the argument DATASET. The microarray platform name PLATFORM
%    must be specified, and must be one of the microarray names found in
%    "ontologies.uarray_platforms".
%
%    IMPORT_UARRAY_DATA(..., 'SampleMap', SAMPLE_MAP) uses the associative map
%    SAMPLE_MAP to specify sample IDs for samples based on their filename. The
%    default is to give no sample IDs for the samples.

% Author: Matti Annala <matti.annala@tut.fi>

function raw = import_uarray_raw(platform, varargin)

files = {};

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'IncludeFiles')
		files = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

metadata.type = 'Microarray probe intensities';
	
import_method = '';

fprintf(1, 'Loading microarray probe definitions...\n');

if strcmpi(platform, 'MDA custom aCGH')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/mda_custom_acgh/probes']);
elseif strcmpi(platform, 'Agilent Human GE')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_human_ge_v1/probes']);
elseif strcmpi(platform, 'Agilent Human GE v2')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_human_ge_v2/probes']);
elseif strcmpi(platform, 'Agilent 244K TCGA custom 1')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_244k_tcga_custom_01/probes']);
elseif strcmpi(platform, 'Agilent 244K TCGA custom 2')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_244k_tcga_custom_02/probes']);
elseif strcmpi(platform, 'Agilent 244K TCGA custom 3')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_244k_tcga_custom_03/probes']);
elseif rx(platform, 'agilent.*244a')       % Agilent HG CGH 244A
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_hg_cgh_244a/probes']);
elseif rx(platform, 'agilent.*cgh.*44b')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_hg_cgh_44b/probes']);
elseif strcmpi(platform, 'Agilent G3 Human CGH 1M')
	import_method = 'agilent_dual';
	probes = load([ppath '/platforms/agilent_g3_cgh_1m/probes']);
elseif strcmpi(platform, 'Agilent Human miRNA 8x15K')
	import_method = 'agilent_single';
	probes = load([ppath '/platforms/agilent_human_mirna_8x15k/probes']);
elseif strcmpi(platform, 'Agilent Human miRNA 8x15K v2')
	import_method = 'agilent_single';
	probes = load([ppath '/platforms/agilent_human_mirna_8x15k_v2/probes']);
elseif strcmpi(platform, 'Agilent Human miRNA 8x15K v3')
	import_method = 'agilent_single';
	probes = load([ppath '/platforms/agilent_human_mirna_8x15k_v3/probes']);
elseif strcmpi(platform, 'Affymetrix HG U133A')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_hg_u133a/probes']);
elseif strcmpi(platform, 'Affymetrix U95A v2')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_u95a_v2/probes']);
elseif strcmpi(platform, 'Affymetrix U95B')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_u95b/probes']);
elseif strcmpi(platform, 'Affymetrix U95C')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_u95c/probes']);
elseif strcmpi(platform, 'Affymetrix HT HG U133A')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_ht_hg_u133a/probes']);
elseif strcmpi(platform, 'Affymetrix HG U133 Plus 2.0')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_hg_u133_plus_2/probes']);
elseif strcmpi(platform, 'Affymetrix Human Exon 1.0 ST')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_huex_1_0_st/probes']);
elseif strcmpi(platform, 'Affymetrix GW SNP 6')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_gw_snp_6/probes']);
elseif strcmpi(platform, 'Affymetrix Human Mapping 250K Sty')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_mapping_250k_sty/probes']);
elseif strcmpi(platform, 'Merja siru')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/merja_siru/probes']);
elseif strcmpi(platform, 'Affymetrix Human Mapping 250K Nsp')
	import_method = 'affy_cel';
	probes = load([ppath '/platforms/affy_mapping_250k_nsp/probes']);
else
	error 'Unrecognized array platform specified.';
end

probes = probes.probes;

uarray_regexp = '^(.+)\.(cel|txt)(\.bz2|\.gz)?$';
compressed_regexp = '^(.+)(\.bz2|\.gz)$';

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
	end
	files = files(:);
end

if length(files) == 0
	fprintf(1, 'No microarray samples found.\n');
	return;
end

S = 0;
raw_mean = {};

for k = 1:length(files)
	
	file = files{k};
	if exist(files{k}) == 0
		fprintf(1, 'File %s does not exist. Skipping it...\n', files{k});
		failed(k) = true;
		continue;
	end
	
	fprintf(1, 'Importing %s...\n', file);
	sample_id = regexprep(file, '\.(cel|txt).*', '', 'ignorecase');
	
	if rx(file, compressed_regexp)
		decompressed_file = decompress_pipe(file);
	else
		decompressed_file = file;
	end
	
	try
		data = struct;
		if strcmp(import_method, 'agilent_dual')
			channels = read_uarray_sample_agilent(decompressed_file, probes, ...
				{'gProcessedSignal', 'rProcessedSignal'});
			data.green_mean = channels(:, 1);
			data.red_mean = channels(:, 2);
		elseif strcmp(import_method, 'agilent_single')
			data.mean = read_uarray_sample_agilent(decompressed_file, probes);
		elseif strcmp(import_method, 'affy_cel')
			data = read_uarray_sample_cel(decompressed_file, probes);
		end
	catch exception
		fprintf(1, 'Skipping sample due to import failure. Reason:\n');
		fprintf(1, '%s\n', exception.message);
		continue;
	end
	
	if strcmp(import_method, 'agilent_dual')
		S = S + 1;
		raw_mean{S} = data.green_mean;
		metadata.sample_id{S} = [sample_id '(Cy3)'];
		metadata.uarray_filename{S} = file;
		metadata.uarray_channel{S} = 'Cy3';
		
		S = S + 1;
		raw_mean{S} = data.red_mean;
		metadata.sample_id{S} = [sample_id '(Cy5)'];
		metadata.uarray_filename{S} = file;
		metadata.uarray_channel{S} = 'Cy5';
	else
		S = S + 1;
		raw_mean{S} = data.mean;
		metadata.sample_id{S} = sample_id;
		metadata.uarray_filename{S} = file;
	end
end

metadata.platform = repmat({platform}, 1, S);

raw.mean = zeros(length(raw_mean{1}), S);
for s = 1:S, raw.mean(:, s) = raw_mean{s}; end
raw.meta = metadata;


