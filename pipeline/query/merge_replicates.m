
% MERGE_REPLICATES   Produce a consensus from technical replicate measurements
%
%    MERGED = MERGE_REPLICATES(DATA) looks for technical replicates within
%    the DATA, and builds consensus samples from the replicates. The way in
%    which the consensus samples are built depends on the dataset type. The
%    function considers two or more samples to be technical replicates if they
%    share the same sample ID. Samples with an unknown ID ('-' or empty) are
%    never merged.
%
%    By default, the function will merge and collapse N technical replicates
%    with the same sample ID into a single sample, taking the sample metadata
%    from the first sample (the assumption is that the metadata is the same
%    for technical replicates).
%
%    MERGE_REPLICATES(..., 'Collapse', false) tells the function to not collapse
%    merged technical replicates into one sample, but to instead replace each
%    replicate sample by the consensus sample. The number of samples and the
%    metadata for each sample are unchanged. Default is to collapse replicates.
%
%    Consensus method for different data types:
%    - Quantile normalization + median across replicates:
%        * Microarray probe intensities
%        * Gene expression
%        * miRNA expression
%    - Union:
%        * Mutations

% Author: Matti Annala <matti.annala@tut.fi>

function merged = merge_replicates(data, varargin)

collapse = true;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Collapse')
		collapse = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

if ~isfield(data, 'Meta') || ~isfield(data.Meta, 'Type') || ...
	~isfield(data.Meta, 'Sample') || ~isfield(data.Meta.Sample, 'ID')
	fprintf(1, 'Dataset has no sample IDs. Cannot merge replicates.\n');
	merged = data;
end

% We must make sure that we don't merge samples with unknown sample IDs. The
% simplest way to ensure this is to give these samples temporary unique IDs.
for s = 1:length(data.Meta.Sample.ID)
	unique_id = 1;
	if strcmp(data.Meta.Sample.ID{s}, '-') || length(data.Meta.Sample.ID{s})==0
		data.Meta.Sample.ID{s} = sprintf('__temp_%d', unique_id);
		unique_id = unique_id + 1;
	end
end

[uniq_samples, ~, groups] = unique(data.Meta.Sample.ID);
replicates = cell(length(uniq_samples), 1);
sample_perm = zeros(length(uniq_samples), 1);

for k = 1:length(uniq_samples)
	rep = find(groups == k);
	replicates{k, 1} = rep;
	sample_perm(k) = rep(1);
end

merged = struct;

if collapse
	merged.Meta = filter_struct(data.Meta, sample_perm);
else
	merged.Meta = data.Meta;
end
	
if strcmp(data.Meta.Type, 'Microarray probe intensities') || ...
	strcmp(data.Meta.Type, 'Gene expression') || ...
	strcmp(data.Meta.Type, 'miRNA expression')
	
	progress = Progress;

	for r = 1:length(replicates)
		rep = replicates{r};
		qnorm = quantilenorm(data.Mean(:, rep));
		merged.Mean = zeros(size(data.Mean, 1), length(uniq_samples));
		merged.Mean(:, r) = median(qnorm, 2);
		progress.update(r / length(replicates));
	end
	
	if ~collapse
		merged.Mean = merged.Mean(:, groups);
	end
	
elseif strcmp(data.Meta.Type, 'Mutations')
	merged.Meta = rmfield(merged.Meta, 'Ref');
	
	for r = 1:length(replicates)
		rep = replicates{r};
		if length(rep) == 1
			merged.Mutations{1, r} = data.Mutations{rep(1)};
		else
			merged.Mutations{1, r} = cat_structs(data.Mutations{rep});
		end
	end
	
	if ~collapse
		merged.Mutations = merged.Mutations(:, groups);
	end
	
else
	error 'Dataset is of an unrecognized type. Cannot merge replicates.';
end


% Remove the temporary unique IDs that were used for samples with unknown
% sample IDs.
for s = 1:length(data.Meta.Sample.ID)
	if length(data.Meta.Sample.ID{s} >= 7) && ...
		strcmp(data.Meta.Sample.ID{s}, '__temp_')
		data.Meta.Sample.ID{s} = '-';
	end
end


