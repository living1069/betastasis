
% Author: Matti Annala <matti.annala@tut.fi>

function alignments = import_alignments(varargin)

%tmp = temporary('alignments');

recursive = true;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Recursive')
		recursive = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

files = find_files('\.bam$');
files = sort_nat(files);

S = length(files);   % Number of sequence samples found.

files = cellfun(@absolutepath, files, 'UniformOutput', false);

alignments.meta.type = 'Sequence alignments';
alignments.meta.sample_id = regexprep(files', '.*/([^/])+.bam', '$1');
alignments.url = strrep(files', '.bam', '');
alignments.format = repmat({'BAM'}, 1, S);

