function [] = seq_error_rate(reads, image)

seq_files = seq_resource_files(reads);

[color, quality] = seq_read_type(seq_files{1});

fprintf(1, 'Aligning reads to transcriptome using Bowtie...\n');
[alignments_tmp, out] = bowtie_align(seq_files{1}, 'transcripts', ...
	'-p4 -v1 -m1 --suppress 1,2,3,4,5,7');

fprintf(1, 'Reading aligned reads into memory...\n');
alignments_file = fopen(alignments_tmp);
data = textscan(alignments_file, '%s %s', 'Delimiter', '\t');
fclose(alignments_file);

%fprintf(1, 'Alignments stored at %s.\n', alignments_tmp);
system(['rm ' alignments_tmp]);
%return;

qualities = data{1};
mismatches = data{2};
clear data;

if length(qualities) ~= length(mismatches)
    error 'Invalid alignment data.';
end

mismatch_count_dist = zeros(1, 10);
read_len_dist = zeros(1, 100);
mismatch_offset_dist = zeros(1, 100);

if color
    for k = 1:length(qualities)
    	bad = find(qualities{k} == '!');
    	if isempty(bad)
    		mismatch_count_dist(1) = mismatch_count_dist(1) + 1;
    		continue;
    	end
    	
    	mismatch_count_dist(2) = mismatch_count_dist(2) + 1;
    	mismatch_offset_dist(bad(1)) = mismatch_offset_dist(bad(1)) + 1;
    end
else
    for k = 1:length(qualities)
        read_len = length(qualities{k});
        read_len_dist(read_len) = read_len_dist(read_len) + 1;
        
        mismatch = mismatches{k};
        bad = str2num(mismatch(1:find(mismatch == ':')));
        
    end
end


mismatch_offset_dist = mismatch_offset_dist / length(qualities);

if nargin == 2
	fprintf(1, 'Rendering the mismatch offset distribution...\n');
	figure; bar(mismatch_offset_dist);
	xlim([0 (max(find(mismatch_offset_dist ~= 0)) + 1)]);
	xlabel('Nucleotide offset'); ylabel('Probability of read error');
	saveas(gcf, image);
end

