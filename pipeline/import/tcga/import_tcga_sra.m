function [] = import_tcga_sra(out_dir, varargin)

delete_originals = false;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'delete')
		delete_originals = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end


if nargin < 1, error 'Please provide an output directory.'; end
	
if ~exist(out_dir), [~, ~] = mkdir(out_dir); end

[~, out] = unix('find -name *.lite.sra');
files = textscan(out, '%s', 'Delimiter', '\n');
files = files{1};

fprintf(1, 'Found %d SRA files.\n', length(files));

for f = 1:length(files)
	tokens = regexp(files{f}, '/(TCGA-[^/-]+-[^/-]+-[^/-]+-[^/]+)/', 'tokens');
	token = tokens{1}; sample_id = token{1};
	
	tokens = regexp(files{f}, '(SRR\d+)', 'tokens');
	token = tokens{1}; srr_id = token{1};

	[status, out] = unix(sprintf('fastq-dump --fasta -O %s %s', ...
		out_dir, files{f}));
	if status ~= 0, error('fastq-dump failed:\n%s\n', out); end
		
	if exist([out_dir '/' srr_id '.fasta'])
		unix(sprintf('cat %s >> %s && rm %s', [out_dir '/' srr_id '.fasta'], ...
			[out_dir '/' sample_id '.fa'], [out_dir '/' srr_id '.fasta']));
	elseif exist([out_dir '/' srr_id '_1.fasta'])
		unix(sprintf('cat %s >> %s && rm %s', ...
			[out_dir '/' srr_id '_1.fasta'], ...
			[out_dir '/' sample_id '_1.fa'], ...
			[out_dir '/' srr_id '_1.fasta']));
		unix(sprintf('cat %s >> %s && rm %s', ...
			[out_dir '/' srr_id '_2.fasta'], ...
			[out_dir '/' sample_id '_2.fa'], ...
			[out_dir '/' srr_id '_2.fasta']));
	end
	
	if delete_originals, delete(files{f}); end
end

