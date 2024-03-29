global pipeline_config;
global pipeline_root;

pipeline_root = '/data/csb/pipeline';

warning off MATLAB:dispatcher:nameConflict;

src_root = strrep(which('pipeline'), [filesep 'pipeline.m'], '');
d = strcat([src_root '/'], { ...
	'config' 'export' 'fmatrix' 'fusion' ...
	'import' 'import/dbsnp' 'import/ensembl' 'import/microarray' 'import/mirbase' ...
	'import/refseq' 'import/tcga' 'import/ucsc' ...
	'microarray' 'microarray/cgh' 'mirna' ...
	'misc' 'ontologies' 'organism' 'query' ...
	'sequencing' 'sequencing/aligners' 'sequencing/cnv_seq' ...
	'sequencing/microbes' 'sequencing/variation' 'sequencing/transform' ...
	'sysbio' ...
});
addpath(d{:});

pipeline_config = struct;
if ~isempty(pipeline_root)
	pipeline_config.RootPath = pipeline_root;
else
	pipeline_config.RootPath = regexprep(src_root, '(/|\\)sources$', '');
end

pipeline_config.MaxThreads = 8;

if isunix || ismac
	pipeline_config.TempDir = ['/data/tmp/' getenv('USER')];
else
	pipeline_config.TempDir = regexprep(tempname, '^(.+)\\(.+?)$', '$1');
end

select_organism('Homo sapiens', '2009');
global organism;

clear pipeline_config src_root d pipeline_root;

if ~isempty(strfind(version, '2010'))
	fprintf('Please hit ENTER to see the command prompt (this is a known Matlab 2010 bug).\n');
end


