
% Author: Matti Annala <matti.annala@tut.fi>

function [] = tophat_fusion(reads)

global organism;

tmp = temporary('tophat_fusion')

max_mismatches = 1;
min_distance = 1e6;
mate_inner_distance = 500;
mate_inner_distance_stdev = 200;
fusion_anchor_len = 20;

index = sprintf('%s/tools/bowtie-indexes/homo_sapiens/2009/genome', ppath);

S = length(reads.url);

if any(rx(reads.paired, 'single'))
	error 'Tophat-fusion analysis currently only supports paired end reads.';
end

flags = '-p8 --no-coverage-search';
flags = sprintf('%s -r %d --mate-std-dev %d', flags, mate_inner_distance, ...
	mate_inner_distance_stdev);
flags = sprintf('%s --fusion-min-dist %d --fusion-read-mismatches %d', ...
	flags, min_distance, max_mismatches);
flags = sprintf('%s --fusion-anchor-length %d', flags, fusion_anchor_len);
flags = sprintf('%s --fusion-multireads 1 --fusion-multipairs 1', flags);


for s = 1:S
	fprintf(1, 'Running tophat-fusion on sample %s...\n', ...
		reads.meta.sample_id{s});
	
	extracted = extract_reads(filter(reads, s), 'UsePipe', false);
	
	out_dir = [tmp reads.meta.sample_id{s}];
	[~, ~] = mkdir(out_dir);
	cwd = cd(out_dir);
	
	[status, out] = unix(sprintf( ...
		'%s/tools/tophat-fusion/tophat-fusion %s %s %s %s', ...
		ppath, flags, index, [extracted.url{1} '_1.fa'], ...
		[extracted.url{1} '_2.fa']));
	if status ~= 0, error('Tophat-fusion returned an error:\n%s\n', out); end
	
	cd(cwd);
end


