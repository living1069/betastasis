
function txome_assembly = assemble_txome_cufflinks(reads, varargin)

global organism;
global pipeline_config;

max_threads = pipeline_config.MaxThreads;

tmp = temporary('tophat')

S = length(reads.url);

if any(~rx(reads.format, 'FASTA'))
	error 'Only FASTA reads are supported by Cufflinks at the moment.';
end

for s = 1:S
	out_dir = [tmp '/' reads.meta.sample_id{s}];
	[~, ~] = mkdir(out_dir);
	
	tophat_flags = '';
	if rx(reads.space{1}, 'color'), tophat_flags = '-C'; end
	tophat_flags = sprintf('%s -p%d', tophat_flags, max_threads);
	
	extracted = extract_reads(filter(reads, s));
	
	if rx(reads.paired, 'single')
		read_files = [reads.url{1} '_1.fa'];
	elseif rx(reads.paired, 'paired')
		read_files = [reads.url{1} '_1.fa ' reads.url{1} '_2.fa'];
	end
	
	status = unix(sprintf('tophat -o %s ' ...
		'%s/tools/bowtie/indexes/%s/%s/genome %s', ...
		out_dir, ppath, organism.Name, organism.Version, read_files));
	if status ~= 0, error 'Tophat returned an error.'; end
	
	
end

txome_assembly = [];
return;



