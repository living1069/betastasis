function coverage = seq_expression(reads, features, varargin)

max_mismatches = 1;

if ischar(features)
	seq = features;
	features = struct;
	features.name = {'Unnamed'};
	features.sequence = {seq};
elseif isfield(features, 'name') && isfield(features, 'sequence')

else
	error(['Sequence features must be provided as a structure with ' ...
		'fields "name" and "sequence".']);
end

S = length(reads.url);
F = length(features.sequence);

coverage.reads = cell(F, S);
coverage.readcou = zeros(F, S);

for s = 1:S
	alignments = bowtie2_align(filter(reads, s), features, ...
		sprintf('-v%d -k10 -m10', max_mismatches));
		
	al = all_alignments(alignments);
	al.target
	al.sequence

	for f = 1:F
		coverage.read_count(f, s) = sum(strcmp(features.name{f}, al.target));
		coverage.reads{f, s} = al.sequence(strcmp(features.name{f}, al.target));
	end
end

