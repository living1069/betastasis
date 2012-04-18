function [alignments, unaligned] = bowtie_align(reads, index, user_flags, ...
	varargin)

global pipeline_config;

tmp = temporary('bowtie_align');

ignore_pairs = true;
max_threads = pipeline_config.MaxThreads;
output = 'Bowtie';

for k = 1:2:length(varargin)
	if regexpi(varargin{k}, 'ignore.*pair')
		ignore_pairs = varargin{k+1}; continue;
	end
	if regexpi(varargin{k}, 'output')
		output = varargin{k+1};
		if output == false, output = 'none'; end
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

S = length(reads.url);

unaligned = [];
if nargout == 2
	unaligned_tmp = temporary('bowtie_unaligned');
	unaligned = reads;
	unaligned.url = strcat(unaligned_tmp, reads.meta.sample_id);
	unaligned.format = repmat({ 'FASTA (gzip)' }, 1, S);
end

alignments = struct;
alignments.meta = reads.meta;
alignments.url = strcat(tmp, reads.meta.sample_id);
alignments.format = repmat({'Bowtie'}, 1, S);
alignments.paired = reads.paired;
if ignore_pairs
	alignments.paired = repmat({'Single'}, 1, S);
end






% We currently require that all read files must be in the same space.
if length(unique(reads.space)) > 1, error 'All reads must be in same space'; end
color = rx(reads.space, 'color');

if ~rx(output, 'none|bowtie')
	error 'Only Bowtie format output is currently supported.';
end





color_option = '';
index_suffix = '';
if color
	color_option = '-C';
	index_suffix = '_colorspace';
end


% Determine the type of index that the user wishes to align against.
if ischar(index)
	% If the index name is specified as a relative path, prefix it with the
	% pipeline directory that holds all Bowtie indices.
	if index(1) ~= '/'
		index = bowtie_index(index);
	end
	index_name = [index index_suffix];
else
	fasta_tmp = [tmp 'index.fasta'];
	index_tmp = [tmp 'index'];

	write_seq_fasta(index, fasta_tmp);

	index_name = [index_tmp index_suffix];

	[status, out] = unix(sprintf('bowtie-build %s %s %s', ...
		color_option, fasta_tmp, index_name));
	if status ~= 0, error('Bowtie index construction failed:\n%s\n.', out);end
end



if any(~rx(reads.format, 'FASTA'))
	error 'Only FASTA files can be aligned at the moment.';
end







for s = 1:S
	
	paired = rx(reads.paired{s}, 'paired');
	index_suffix = '';
	flags = [user_flags ' -f -B1'];

	% Determine general alignment flags.
	if color
		flags = [flags ' -C'];
		index_suffix = '_colorspace';
	end

	flags = sprintf('%s -p%d', flags, max_threads);

	if rx(output, 'none')
		flags = [flags ' --suppress 1,2,3,4,5,6,7,8'];
	else
		flags = [flags ' --suppress 6,7,8'];
	end

	fprintf(1, '-> Invoking Bowtie with flags "%s"...\n', flags);

	if paired && ignore_pairs == false
		% The user wishes to perform proper paired end alignment.
		error 'Proper paired end alignment not supported yet.';

	elseif paired && ignore_pairs
		% FIXME: We assume that reads have /1 and /2 suffixes.
		
		if ~isempty(unaligned)
			% Write unaligned reads into a pipe, run separate_fasta_pairs.py
			% on data coming out of the pipe, and feed the separated pairs
			% into two more pipes that will compress the data on the fly.
			unaligned_tmp = [ptemp '.pair_split_pipe'];
			unix(sprintf( ...
				['mkfifo %s && python %s/sources/sequencing/transform/' ...
				'separate_fasta_pairs.py %s %s %s &'], unaligned_tmp, ppath, ...
				unaligned_tmp, compress_pipe([unaligned.url{s} '_1.fa.gz']), ...
				compress_pipe([unaligned.url{s} '_2.fa.gz'])));

			flags = sprintf('%s --un %s', flags, unaligned_tmp);
		end
		
		inputs = strcat(reads.url{s}, {'_1.fa', '_2.fa'});
		if rx(reads.format{s}, 'gzip')
			for k = 1:length(inputs)
				inputs{k} = sprintf('<(gunzip -c %s.gz)', inputs{k});
			end
		end
		
		[status, out] = unix(sprintf('cat %s %s | bowtie %s %s - > %s', ...
			inputs{1}, inputs{2}, flags, index_name, ...
			[alignments.url{s} '.bowtie']));
		if status ~= 0, error('Bowtie read alignment failed:\n%s', out); end
			
		[alignments.total_reads(s), alignments.aligned_reads(s), ...
			alignments.total_alignments(s)] = parse_bowtie_stats(out);
		
	else
		if ~isempty(unaligned)
			flags = sprintf('%s --un %s', flags, ...
				compress_pipe([unaligned.url{s} '.fa.gz']));
		end
		
		input = [reads.url{s} '.fa'];
		if rx(reads.format{s}, 'gzip')
			input = sprintf('<(gunzip -c %s.gz)', input);
		end

		[status, out] = unix(sprintf('bowtie %s %s %s > %s', ...
			flags, index_name, input, [alignments.url{s} '.bowtie']));
		if status ~= 0, error('Bowtie read alignment failed:\n%s', out); end
		
		[alignments.total_reads(s), alignments.aligned_reads(s), ...
			alignments.total_alignments(s)] = parse_bowtie_stats(out);
	end
end







function [total, aligned, alignments] = parse_bowtie_stats(out)

total = NaN;
aligned = NaN;
alignments = NaN;

out = strread(out, '%s', 'delimiter', '\n');
for k = 1:length(out)
	match = regexp(out{k}, 'reads processed: (\d+)', 'tokens');
	if length(match) == 1
		tmp = match{1}; total = str2double(tmp{1});
		continue;
	end
	
	match = regexp(out{k}, 'reads with .* reported alignment: (\d+)', 'tokens');
	if length(match) == 1
		tmp = match{1}; aligned = str2double(tmp{1});
		continue;
	end
	
	match = regexp(out{k}, '^Reported (\d+) alignments to', 'tokens');
	if length(match) == 1
		tmp = match{1}; alignments = str2double(tmp{1});
		continue;
	end
end

if ~isnan(total) && ~isnan(aligned)
	fprintf(1, '-> %d / %d (%.1f%%) reads with at least one alignment\n', ...
		aligned, total, aligned / total * 100);
end


