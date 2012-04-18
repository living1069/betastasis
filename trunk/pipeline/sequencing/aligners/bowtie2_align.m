function [alignments, unaligned] = bowtie2_align(reads, index, flags, varargin)

global organism;
global pipeline_config;

tmp = temporary('bowtie2_align');
unaligned_tmp = temporary('bowtie2_unaligned');

max_threads = pipeline_config.MaxThreads;
ignore_pairs = true;
output = 'BAM';

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'ignore.*pair')
		ignore_pairs = varargin{k+1}; continue;
	end
	if rx(varargin{k}, 'output')
		output = varargin{k+1}; continue;
	end
	if rx(varargin{k}, 'max.*thread')
		max_threads = varargin{k+1}; continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

S = length(reads.url);

unaligned = [];
if nargout == 2
	unaligned = reads;
	unaligned.url = strcat(unaligned_tmp, reads.meta.sample_id);
	unaligned.format = repmat({ 'FASTA (gzip)' }, 1, S);
	unaligned.paired = repmat({ 'Single' }, 1, S);
end

alignments = struct;
alignments.meta = reads.meta;
alignments.url = strcat(tmp, reads.meta.sample_id);
alignments.format = repmat({'BAM'}, 1, S);
alignments.paired = reads.paired;
if ignore_pairs
	alignments.paired = repmat({'Single'}, 1, S);
end






% We currently require that all read files must be in the same space.
if any(rx(reads.space, 'color'))
	error 'Colorspace reads are not supported by Bowtie2.';
end
if any(~rx(reads.format, 'FAST[AQ]'))
	error 'Only FASTA/FASTQ files can be aligned at the moment.';
end
if ~rx(output, 'BAM')
	error 'Only BAM output is currently supported.';
end






% Determine the type of index that the user wishes to align against.
if ischar(index)
	% If the index name is specified as a relative path, prefix it with the
	% pipeline directory that holds all Bowtie indices.
	if index(1) ~= '/'
		index = lower(index);
		if ~rx(index, '^(genome|transcripts|exons|mirnas|pre_mirnas)$')
			error('Bowtie2 index "%s" is not supported.', index);
		end

		index = ['/home/csbgroup/tools/bowtie2-indexes/' ...
			flatten_str(organism.Name) '/' ...
			flatten_str(organism.Version) '/' index];
	end
	index_name = index;
else
	fasta_tmp = [tmp 'index.fasta'];
	index_name = [tmp 'index'];

	write_seq_fasta(index, fasta_tmp);

	[status, out] = unix(sprintf('bowtie2-build %s %s', fasta_tmp, index_name));
	if status ~= 0, error('Bowtie2 index construction failed:\n%s\n.', out);end
end










user_flags = flags;

for s = 1:S
	read_paths = access_reads(filter(reads, s));
	paired = rx(reads.paired{s}, 'paired');
	
	if rx(reads.format{s}, 'FASTQ')
		flags = '';
	elseif rx(reads.format{s}, 'FASTA')
		flags = '-f';
	end
	
	if ~rx(user_flags, '-p\d+')
		flags = sprintf('%s -p%d', flags, max_threads);
	end

	flags = [flags ' ' user_flags];
	
	fprintf('-> Invoking Bowtie2 with flags "%s"...\n', flags);
	
	output_file = [alignments.url{s} '.bam'];
	if rx(output, 'none'), output_file = '/dev/null'; end

	if paired && ignore_pairs == false
		% The user wishes to perform proper paired end alignment.
		error 'Proper paired end alignment not supported yet.';

	elseif paired && ignore_pairs == true
		% THIS IS THE CASE WHERE WE ALIGN LEFT AND RIGHT READS SEPARATELY
		% AND THROW ALL THE ALIGNMENTS INTO ONE SINGLE FILE.
		% FIXME: We assume that reads have /1 and /2 suffixes.
		
		if ~isempty(unaligned)
			flags = sprintf('%s --un %s', flags, ...
				compress_pipe([unaligned.url{s} '.fa.gz']));
		end
		
		[status, out] = unix(sprintf( ...
			['cat %s %s | bowtie2 %s -x %s - | ' ...
			'%s/tools/samtools/samtools view -Sb -o %s -'], ...
			read_paths{1}, read_paths{2}, ...
			flags, index_name, ppath, output_file));
		if status ~= 0, error('Bowtie2 read alignment failed:\n%s\n', out); end
			
		[alignments.total_reads(s), alignments.aligned_reads(s), ...
			alignments.total_alignments(s)] = parse_bowtie_stats(out);
		
	else
		% Single end read alignment.
		if ~isempty(unaligned)
			flags = sprintf('%s --un %s', flags, [unaligned.url{s} '.fa']);
		end

		[status, out] = unix(sprintf( ...
			['bowtie2 %s -x %s %s | ' ...
			'%s/tools/samtools/samtools view -Sb -o %s -'], ...
			flags, index_name, read_paths{1}, ppath, output_file));
		if status ~= 0, error('Bowtie2 read alignment failed:\n%s\n', out); end
		
		[alignments.total_reads(s), alignments.aligned_reads(s), ...
			alignments.total_alignments(s)] = parse_bowtie_stats(out);
	end
end







function [total, aligned, alignments] = parse_bowtie_stats(out)

total = NaN;
aligned = 0;
alignments = 0;

out = strread(out, '%s', 'delimiter', '\n');
for k = 1:length(out)
	match = regexp(out{k}, '(\d+) reads; of these:', 'tokens');
	if length(match) == 1
		tmp = match{1}; total = str2double(tmp{1});
		continue;
	end
	
	match = regexp(out{k}, '(\d+) .* aligned 0 times', 'tokens');
	if length(match) == 1
		tmp = match{1}; aligned = total - str2double(tmp{1});
		continue;
	end
end

if ~isnan(total) && ~isnan(aligned)
	fprintf(1, '-> %d / %d (%.1f%%) reads with at least one alignment\n', ...
		aligned, total, aligned / total * 100);
end


