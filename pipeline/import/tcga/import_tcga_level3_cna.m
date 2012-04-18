
function segments = import_tcga_level3_cna(files)

global organism;

if nargin < 1
	files = find_files('.*_copy_number_analysis.txt');
end

segments = struct;
segments.meta.type = 'Copy number segments';
segments.meta.sample_id = {};
segments.segments = {};

for f = 1:length(files)
	fid = fopen(files{f});
	[data, headers] = readtable(fid, 'Numeric', {'start', 'stop', 'seg.mean'});
	fclose(fid);

	samples = data{1};
	chromosomes = chromosome_sym2num(data{2});
	seg_pos = [data{3} data{4}];
	seg_mean = data{6};
	
	[~, sample_starts] = unique(samples, 'first');
	sample_starts = [sample_starts; length(samples)+1];
	
	for s = 1:length(sample_starts)-1
		segments.meta.sample_id{end+1} = samples{s};
		S = length(segments.meta.sample_id);
		
		segments.genome_build{S} = '';
		
		lines = sample_starts(s):sample_starts(s+1)-1;
		
		segs = struct;
		segs.chromosome = chromosomes(lines);
		segs.position = seg_pos(lines, :);
		segs.logratio = seg_mean(lines);
		
		segments.segments{S} = segs;
	end
end

