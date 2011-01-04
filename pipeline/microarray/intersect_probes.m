
% INTERSECT_PROBES    Intersect multiple microarray probe definitions
%
%    [INTERSECTED, INCLUDED_PROBES] = INTERSECT_PROBES(ARRAYS) takes the 
%    intersection among the microarray probe definitions in ARRAYS, and
%    constructs new probe definitions that are returned in INTERSECTED.
%    ARRAYS must be a cell array of microarray probe definition structures, 
%    and INTERSECTED is also returned as a cell array of the same size.
%    The function also returns the matrix INCLUDED_PROBES, the columns of
%    which indicate the probes that were included in the intersected probe
%    definition from each array.

% Author: Matti Annala <matti.annala@tut.fi>

function [intersected, included_probes] = intersect_probes(arrays, varargin)

check_position = false;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'CheckPosition')
		check_position = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

for a = 1:length(arrays)
	arrays{a}.Sequence = upper(arrays{a}.Sequence);
end

common_seq = arrays{1}.Sequence;
for a = 2:length(arrays)
	common_seq = intersect(common_seq, arrays{a}.Sequence);
end

array_sorted_probes = cell(length(arrays), 1);
array_idx = cell(length(arrays), 1);
for a = 1:length(arrays)
	[array_sorted_probes{a}, array_idx{a}] = sort(arrays{a}.Sequence);
end

progress = Progress;

S = length(common_seq);
pos = ones(1, length(arrays));

% For each probe sequence that is present on at least one of the arrays, check
% if that probe sequence is present on all of the arrays, and how many common
% copies the arrays have. For instance, if array A has 3 copies of the probe
% sequence and array B has 2 copies, then the intersected list of probes will
% contain 2 copies.

included_probes = zeros(length(arrays{1}.Sequence), length(arrays));
ipos = 1;

for s = 1:S
	seq = common_seq{s};
	
	num_seq_probes = zeros(1, length(arrays));
	
	for a = 1:length(arrays)
		while ~strcmpi(array_sorted_probes{a}{pos(a)}, seq)
			pos(a) = pos(a) + 1;
		end
		
		while pos(a) <= length(array_sorted_probes{a}) && ...
			strcmpi(array_sorted_probes{a}{pos(a)}, seq)
			
			num_seq_probes(a) = num_seq_probes(a) + 1;
			pos(a) = pos(a) + 1;
		end
	end
	
	min_copies = min(num_seq_probes);
	
	for a = 1:length(arrays)
		included_probes(ipos:ipos+min_copies-1, a) = ...
			array_idx{a}(pos(a)-min_copies:pos(a)-1);
	end
	
	ipos = ipos + min_copies;
	progress.update(s / S);
end

included_probes = included_probes(1:ipos-1, :);

fprintf(1, 'Found %d common probes.\n', size(included_probes, 1));

for a = 1:length(arrays)
	iarray = struct;
	iarray.XPos = arrays{a}.XPos(included_probes(:, a));
	iarray.YPos = arrays{a}.YPos(included_probes(:, a));
	iarray.Sequence = arrays{a}.Sequence(included_probes(:, a));
	
	intersected{a} = iarray;
end

