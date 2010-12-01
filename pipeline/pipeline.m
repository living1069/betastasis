global pipeline_config;
global pipeline_root;

src_root = strrep(which('pipeline'), [filesep 'pipeline.m'], '');
d = strcat([src_root '/'], { ...
	'config' 'fusion' ...
	'import' 'import/dbsnp' 'import/microarray' 'import/refseq' ...
	'import/tcga' 'import/ucsc' ...
	'microarray' 'microarray/cgh' 'mirna' 'misc' ...
	'ontologies' 'organism' 'query' 'sequencing' ...
	'sequencing/aligners' 'sequencing/bowtie' });
addpath(d{:});

pipeline_config = struct;
if ~isempty(pipeline_root)
	pipeline_config.RootPath = pipeline_root;
else
	pipeline_config.RootPath = regexprep(src_root, '(/|\\)sources$', '');
end

pipeline_config.MaxThreads = 4;

pipeline_config.TCGA = struct;
pipeline_config.TCGA.Path = '/worktmp/TCGA';

if isunix || ismac
	pipeline_config.TempDir = '/worktmp/tmp';
else
	pipeline_config.TempDir = regexprep(tempname, '^(.+)\\(.+?)$', '$1');
end

[~, ~] = mkdir([ppath '/datasets']);
pipeline_config.Repositories = { ...
	LocalRepository('local', [ppath '/datasets']) ...
};

select_organism('Homo sapiens', '2009');
global organism;

clear pipeline_config src_root d pipeline_root;

if ~isempty(strfind(version, '2010'))
	fprintf(1, 'Please hit ENTER to see the command prompt (this is a known Matlab 2010 bug).\n');
end

