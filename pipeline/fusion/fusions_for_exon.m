function rearrangements = fusions_for_exon(reads, exon, n_prime, min_read_len)

global organism;
transcripts = organism.Transcripts;
exons = organism.Exons;

if nargin < 3
	min_read_len = NaN;
end

seq_files = seq_resource_files(reads);

[color, quality] = seq_read_type(seq_files{1});

rearrangements = struct;
rearrangements.Fusions = cell(1, length(seq_files));

if n_prime == 5
	prime_str = '5''';
elseif n_prime == 3
	prime_str = '3''';
else
	error 'Invalid specification of exon half within fusion.';
end

fprintf(1, ['Constructing a list of all possible exon-exon junctions with ' ...
            'the given exon at the %s end...\n'], prime_str);

fasta_fname = ptemp();
index_fname = ptemp();

fasta_file = fopen(fasta_fname, 'W');

num_exon_pairs = 0;

fixed_gene = organism.Transcripts.Gene(organism.Exons.Transcript(exon));

for g = [1:fixed_gene-1 fixed_gene+1:length(organism.Genes.Name)]
	
	lens = [];
	ts = organism.Genes.Transcripts(g, 1:organism.Genes.TranscriptCount(g));
	for t = 1:length(ts)
		lens(t) = length(organism.Transcripts.Sequence{ts(t)});
	end
	[~, longest_ts] = max(lens);
	ts = ts(longest_ts(1));
	
	exons = find(organism.Exons.Transcript == ts);

	for ex = 1:length(exons)
		left_exon = left_exons(l_ex);
					
					for r_ex = 1:length(right_exons)
						right_exon = right_exons(r_ex);
						
						left_ts_seq = transcripts.Sequence{ ...
							exons.Transcript(left_exon)};
						right_ts_seq = transcripts.Sequence{ ...
							exons.Transcript(right_exon)};

						left_exon_seq = left_ts_seq( ...
							exons.Position(left_exon, 1): ...
							exons.Position(left_exon, 2));
						right_exon_seq = right_ts_seq( ...
							exons.Position(right_exon, 1): ...
							exons.Position(right_exon, 2));

						if length(left_exon_seq) > 100
							left_exon_seq = left_exon_seq(end-99:end);
						end
						if length(right_exon_seq) > 100
							right_exon_seq = right_exon_seq(1:100);
						end
						
						fprintf(fasta_file, '>%d,%d\n', left_exon, right_exon);
						fprintf(fasta_file, '%s%s\n', ...
							left_exon_seq, right_exon_seq);
							
						num_exon_pairs = num_exon_pairs + 1;
					end
				end
			end
		end
	end
end

fclose(fasta_file);

fprintf(1, 'Number of potential fusion junctions to check: %d\n', ...
	num_exon_pairs);

fprintf(1, ['Constructing a Bowtie index with all exon-exon junctions ' ...
            'between genes of interest...\n']);
	
color_option = '';
index_suffix = '';
if color
	color_option = '-C';
	index_suffix = '_colorspace';
end

[status, ~] = unix(sprintf('%s/tools/bowtie/bowtie-build %s %s %s%s', ...
	ppath, color_option, fasta_fname, index_fname, index_suffix));

if status ~= 0, error 'Bowtie index construction failed.'; end

unaligned_read_files = cell(length(seq_files), 1);
for k = 1:length(seq_files)
	len_filtered_reads_tmp = ptemp();
	unaligned_read_files{k} = ptemp();
	
	% FIXME: Length filtering doesn't currently work with colorspace reads.
	if ~isnan(min_read_len)
		fprintf(1, 'Filtering out reads shorter than %d bp...\n', min_read_len);
		[status, ~] = unix(sprintf( ...
			'%s/sources/sequencing/filter_reads_by_len.py %s %d > %s', ...
			ppath, seq_files{k}, min_read_len, len_filtered_reads_tmp));
	else
		len_filtered_reads_tmp = seq_files{k};
	end
	
	[flags, index_suffix] = bowtie_flags_for_reads(seq_files{k});
	
	fprintf(1, 'Finding unaligned transcriptome reads using Bowtie...\n');
	[status, out] = unix(sprintf(['%s/tools/bowtie/bowtie %s ' ...
		'-p4 -v2 --suppress 5,6,7,8 %s%s %s --un %s > /dev/null'], ...
		ppath, flags, bowtie_index('transcripts'), index_suffix, ...
		len_filtered_reads_tmp, unaligned_read_files{k}));

	if ~isnan(min_read_len)
		system(['rm ' len_filtered_reads_tmp]);
	end
		
	fprintf(1, '%s', out);    % Only stderr is output.
	if status ~= 0, error 'Read alignment failed.'; end
end

for k = 1:length(seq_files)
	rearrangements.Fusions{k} = validate_junctions(unaligned_read_files{k}, ...
		index_fname);
	system(['rm ' unaligned_read_files{k}]);
end

system(['rm ' index_fname '*']);

if isfield(reads, 'Meta')
	rearrangements.Meta = reads.Meta;
else
	rearrangements.Meta = struct();
end

rearrangements.Meta.Type = 'Genetic rearrangements';
rearrangements.Meta.Organism = organism.Name;
rearrangements.Meta.GenomeVersion = organism.GenomeVersion;
	
return;







function validated = validate_junctions(unaligned_reads, index_fname)

global organism;
exons = organism.Exons;

[color, ~] = seq_read_type(unaligned_reads);

fprintf(1, 'Aligning full reads against candidate junctions...\n');
[alignments_tmp, out] = bowtie_align(unaligned_reads, index_fname, ...
	'-p4 -m4 -v2 --suppress 2,6,7,8');

fprintf(1, '%s', out);

fprintf(1, 'Reading aligned reads into memory...\n');
alignments_file = fopen(alignments_tmp);
data = textscan(alignments_file, '%s %s %d %s', 'Delimiter', '\t');
fclose(alignments_file);

system(['rm ' alignments_tmp]);

read_ids = data{1};
junction_exons = data{2};
read_offsets = data{3} + 1;
read_sequences = data{4};
clear data;

% Calculate the lengths of the reference sequences on the 5' side of junctions.
left_lengths = zeros(length(junction_exons), 1);
for k = 1:length(junction_exons)
	ex = sscanf(junction_exons{k}, '%d,%d');
	len = exons.Position(ex(1), 2) - exons.Position(ex(1), 1) + 1;
	if len > 100, len = 100; end
	left_lengths(k) = len;
end

% Only reads that have at least 10 nucleotides on both sides of the junction
% count towards validation.
valid = false(length(read_ids), 1);
for k = 1:length(valid)
	len = length(read_sequences{k});
	offset = read_offsets(k);
	if offset <= left_lengths(k) - 7 && offset + len > left_lengths(k) + 8
		valid(k) = true;
	end
end

read_ids = read_ids(valid);
junction_exons = junction_exons(valid);
read_offsets = read_offsets(valid);
read_sequences = read_sequences(valid);

unique_junctions = unique(junction_exons);

junction_reads_count = zeros(length(unique_junctions), 1);
for k = 1:length(unique_junctions)
	junction_reads_count(k) = sum(strcmp(unique_junctions{k}, junction_exons));
end

[~, sort_indices] = sort(junction_reads_count, 1, 'descend');

validated = struct('Exons', zeros(length(unique_junctions), 2), ...
	'ReadCount', junction_reads_count);

validated.ReadSequences = cell(length(unique_junctions), ...
	max(junction_reads_count));

for k = 1:length(unique_junctions)
	validated.Exons(k, :) = sscanf(unique_junctions{k}, '%d,%d');
	
	indices = find(strcmp(unique_junctions{k}, junction_exons));
	for s = 1:length(indices)
		validated.ReadSequences{k, s} = read_sequences{indices(s)};
	end
end	

return;

