function coverage = find_seq_features(reads, features, varargin)

if ~isfield(features, 'Name') || ~isfield(features, 'Sequence')
	error(['Sequence features must be provided as a structure with fields ' ...
		'"Name" and "Sequence".']);
end

background = [];

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Background')
		background = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);




seq_files = seq_resource_files(reads);
S = length(seq_files);
F = length(features.Sequence);

coverage.Reads = cell(F, S);
coverage.ReadCount = zeros(F, S);

for s = 1:S
	al = align_reads(seq_files{s}, features, ...
		'MaxMismatches', 2, 'AllowAlignments', 10, varargin{:}, ...
		'Columns', 'target,offset,sequence');

	for f = 1:F
		coverage.ReadCount(f, s) = sum(strcmp(features.Name{f}, al.Target));
		coverage.Reads{f, s} = al.Sequence(strcmp(features.Name{f}, al.Target));
	end
end

