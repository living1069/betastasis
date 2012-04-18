function extracted = extract_reads(reads, varargin)

S = length(reads.url);

use_pipes = true;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'pipe')
		use_pipes = varargin{k+1}; continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

tmp = temporary('extract_reads');
extracted = reads;

read_files = seq_filenames(reads);

for s = 1:S
	if rx(reads.format{s}, 'gzip')
		extracted.format{s} = strrep(extracted.format{s}, ' (gzip)', '');
		extracted.url{s} = [tmp reads.meta.sample_id{s}];
		pipe_files = seq_filenames(filter(extracted, s));
		pipe_files = pipe_files{1};
		for f = 1:length(read_files{s})
			if use_pipes
				unix(['mkfifo ' pipe_files{f}]);
				unix(sprintf('gunzip -c %s > %s &', ...
					read_files{s}{f}, pipe_files{f}));
			else
				unix(sprintf('gunzip -c %s > %s', ...
					read_files{s}{f}, pipe_files{f}));
			end
		end
	end
end


