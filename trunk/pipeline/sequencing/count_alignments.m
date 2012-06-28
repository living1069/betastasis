
function expr = count_alignments(alignments, bed_file)

global organism;

num_parallel = 8;

tmp = temporary('count_alignments');

S = length(alignments.url);

% Count number of columns in the BED file.
data = readtable(bed_file, 'Header', false, 'NumLines', 1);
num_annot_cols = length(data);

expr = struct;
expr.meta = alignments.meta;

unix(sprintf('cd %s && split -l 10000 %s', tmp, bed_file));
files = dir(tmp); files = {files.name};
files(~rx(files, '^x')) = [];
files = sort(files);
files = strcat(tmp, files);

bam_urls = '';
for s = 1:S, bam_urls = [bam_urls alignments.url{s} '.bam ']; end
	
job_status = zeros(1, length(files));
while ~all(job_status == 2)
	for f = 1:length(files)
		
		if exist([files{f} '.count'])
			job_status(f) = 2;
		
		elseif job_status(f) == 0 && sum(job_status == 1) < num_parallel
			fprintf('Counting alignments for feature batch %s...\n', ...
				basename(files{f}));
			unix(sprintf([ ...
				'bedtools multicov -bams %s -bed %s > %s.incomplete.count &&'...
				'mv %s.incomplete.count %s.count &'], ...
				bam_urls, files{f}, files{f}, files{f}, files{f}));
			
			job_status(f) = 1;
		end
	end

	pause(10);
end

unix(sprintf('cat %s > %s/total.count', sprintf('%s.count ', files{:}), tmp));

data = readtable([tmp 'total.count'], 'Header', false, ...
	'Numeric', num_annot_cols+1:num_annot_cols+S);
expr.rows.name = data{4};
expr.rows.chromosome = chromosome_sym2num(data{1});
expr.rows.position = [str2double(data{2}), str2double(data{3})];

expr.mean = nan(length(data{1}), S);
for s = 1:S, expr.mean(:, s) = data{num_annot_cols+s}; end

