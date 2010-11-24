
% MERGE_REPLICATES   Produce a consensus from technical replicate measurements
%
%    MERGED = MERGE_REPLICATES(DATA) looks for technical replicates within
%    the DATA, and builds consensus samples from the replicates. The way in
%    which the consensus samples are built depends on the dataset type. The
%    function considers two or more samples to be technical replicates if they
%    share the same sample ID.
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

uniq_samples = unique(data.Meta.Sample.ID);
replicates = cell(length(uniq_samples), 1);
sample_perm = zeros(length(uniq_samples), 1);

for k = 1:length(uniq_samples)
	rep = find(strcmp(uniq_samples{k}, data.Meta.Sample.ID));
	replicates{k, 1} = rep;
	sample_perm(k) = rep(1);
end

merged = struct;

if collapse
	merged.Meta = filter_struct(data.Meta, sample_perm);
	
	if strcmp(data.Meta.Type, 'Microarray probe intensities') || ...
		strcmp(data.Meta.Type, 'Gene expression') || ...
		strcmp(data.Meta.Type, 'miRNA expression')
		
		progress = Progress;
	
		merged.Mean = zeros(size(data.Mean, 1), length(uniq_samples));
		for r = 1:length(replicates)
			rep = replicates{r};
			qnorm = quantilenorm(data.Mean(:, rep));
			merged.Mean(:, r) = median(qnorm, 2);
			progress.update(r / length(replicates));
		end
	else
		error 'Dataset is of an unrecognized type. Cannot merge replicates.';
	end
else
	merged.Meta = data.Meta;
	
	if strcmp(data.Meta.Type, 'Microarray probe intensities') || ...
		strcmp(data.Meta.Type, 'Gene expression') || ...
		strcmp(data.Meta.Type, 'miRNA expression')
	
		progress = Progress;

		merged.Mean = zeros(size(data.Mean, 1), length(uniq_samples));
		for r = 1:length(replicates)
			rep = replicates{r};
			qnorm = quantilenorm(data.Mean(:, rep));
			merged.Mean(:, rep) = repmat(median(qnorm, 2), 1, length(rep));
			progress.update(r / length(replicates));
		end
	else
		error 'Dataset is of an unrecognized type. Cannot merge replicates.';
	end
end


