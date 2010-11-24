function files = seq_resource_files(reads)

% Figure out what kind of format the user specified the sequence files in.
files = {};
if ischar(reads)
	files = { reads };
elseif iscell(reads)
	files = reads;
elseif isstruct(reads) && isfield(reads, 'SequenceResource')
	files = strcat([ppath '/datasets/'], reads.SequenceResource);
else
	error 'Sequence file must be specified as a string, cell array or struct.';
end

