
% PAIRED_SAMPLES    Constructs samples pairs based on two input query sets
%
%    [PA, PB] = PAIRED_SAMPLES(A, B, KEY) takes two query sets as input and
%    constructs a list of sample pairs according to a user defined criterion
%    specified in the argument KEY. The function returns two new query sets
%    PA and PB (of equal size) where the first sample in PA matches with the
%    first sample in PB, and so forth for all samples.
%    
%    This function is useful in experiments where multiple test samples need to
%    be compared with their respective reference samples in a pairwise fashion.
%
%    Supported values for argument KEY:
%    - 'Patient': Pair samples based on patient ID
%    - 'Sample': Pair samples based on sample ID
%    - 'Filename': Pair samples based on original filename
%                  (useful for 2-channel microarrays)
%
%    Example:
%    We have a query set A of three tumor samples from patients A23, A54 and A12
%    (in that order). We also have a query set B of two adjacent normal samples
%    from patients A12 and A54. We then run
%    
%        [PA, PB] = PAIRED_SAMPLES(A, B, 'Patient')
%    
%    to build two new query sets PA and PB where the samples are paired
%    according to the patient ID. The patient IDs for samples in the new query
%    set PA are now A54, A12 (in that order). And for the new query set PB,
%    the patient IDs are also A54, A12, since we paired the samples. Note that
%    the tumor sample for patient A23 was discarded altogether, since no
%    matching reference sample was found in query set B.

% Author: Matti Annala <matti.annala@tut.fi>

function [A, B] = paired_samples(A, B, key)

if nargin > 2
	if regexpi(key, 'patient.?id|patient') 
		key = 'Patient.ID';
	elseif regexpi(key, 'sample.?id|sample')
		key = 'Sample.ID';
	elseif regexpi(key, '(sample.)?filename')
		key = 'Sample.Filename';
	else
		error 'Specified key is not supported.';
	end
else
	key = 'Patient.ID';
end

if isfield(A, 'Meta') ~= isfield(B, 'Meta')
	error 'A and B must either both be query sets or both be realized.';
end

meta_A = A;
meta_B = B;
if isfield(A, 'Meta')
	meta_A = A.Meta;
	meta_B = B.Meta;
end

eval(['keys_a = meta_A.' key ';']);
eval(['keys_b = meta_B.' key ';']);

paired_a = true(length(keys_a), 1);
paired_b = [];

num_ambiguous = 0;

for k = 1:length(keys_a)
	if any(find(strcmp(keys_a{k}, keys_a)) < k)
		paired_a(k) = false;
		continue;
	end
	
	matches = find(strcmp(keys_a{k}, keys_b));
	if length(matches) ~= 1
		if length(matches) > 1
			num_ambiguous = num_ambiguous + 1;
		end
		paired_a(k) = false;
		continue;
	end
	
	paired_b(end + 1) = matches(1);
end

if num_ambiguous > 0
	fprintf(1, ['WARNING: %d samples were discarded because they had more ' ...
		'than one matching pair.\n'], num_ambiguous);
end

if isempty(paired_b)
	fprintf(1, 'WARNING: No paired samples found.\n');
	A = [];
	B = [];
	return;
end

if isfield(A, 'Meta')
	A = filter_realized(A, paired_a);
	B = filter_realized(B, paired_b);
else
	A = filter_query(A, paired_a);
	B = filter_query(B, paired_b);
end






function ds = filter_realized(ds, selected)

fields = fieldnames(ds);
for k = 1:length(fields)
	if strcmp(fields{k}, 'Meta'), continue, end
	
	f = getfield(ds, fields{k});
	if isstruct(f), continue, end
	
	ds = setfield(ds, fields{k}, f(:, selected));
end

ds.Meta = filter_query(ds.Meta, selected);

