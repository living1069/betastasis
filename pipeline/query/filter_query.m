
% FILTER_QUERY         Filter a query set using a metadata predicate.
%
%    QSET = FILTER_QUERY(QSET, PRED) filters the query set QSET so that only
%    samples fulfilling the predicate PRED are left in the query set. The
%    predicate PRED can be specified either as a string in the pipeline's query
%    language, or as a logical or index vector. The latter usage allows you
%    to use normal Matlab expressions and Boolean operators.
%
%    List of valid predicate keys:
%    - SAMPLE ID (numeric)
%    - SAMPLE TYPE (string)
%    - ORGAN / TISSUE (string)
%    - PATIENT ID (string)
%    - GENDER (string)
%    - STATUS (string)
%    - ETHNICITY (string)
%    - RACE (string)
%    - TREATMENT (string)
%    - SURVIVAL TIME (numeric, in days)
%    - CENSORED (0, 1 or NaN)
%    - AGE AT LAST FOLLOWUP (numeric, in years)
%    - AGE AT DEATH (numeric, in years)
%    - AGE AT DIAGNOSIS (numeric, in years)
%    - AGE AT BIRTH (numeric, in years)
%    - FILENAME (string)
%    - CHANNEL (string)
%
%    Valid predicates for metadata of type string:
%    - case-insensitive full string match:
%          'patient = TCGA-02-1032'    OR (equiv.)   'patient == TCGA-02-1032'
%          'patient ~= TCGA-02-1032'   OR (equiv.)   'patient != TCGA-02-1032'
%    - regular expression string match:
%          'sample type ~ tumor'
%          'sample type !~ tumor'
%
%    See also QUERY, REALIZE.

% Author: Matti Annala <matti.annala@tut.fi>

function filtered = filter_query(qset, predicate)

filtered = qset;

if isnumeric(predicate) || islogical(predicate)
	% Assume that we're dealing with an index vector or logical vector.
	filtered = filter_struct(filtered, predicate);
	return;
end

if ~isfield(qset, 'Resource')
	error 'filter_query() only works on query sets.';
end

if isempty(regexp(predicate, '\S')), return, end
	
	
	
	
	
filtered = filter_by_string(qset, 'Patient.ID', 'patient\s*id|patient', ...
	predicate);
if isstruct(filtered), return, end

filtered = filter_by_string(qset, 'Patient.Gender', 'gender|sex', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_string(qset, 'Patient.Status', ...
	'patient\s*status|status|health|state|patient\s*state', predicate);
if isstruct(filtered), return, end

filtered = filter_by_string(qset, 'Patient.Ethnicity', ...
	'ethnicity|nationality|ethnic\s*group', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_string(qset, 'Patient.Race', 'race|color', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_string(qset, 'Patient.Treatment', ...
	['treatments?|drugs?|treatment\s*history|medical\s*history'], predicate);
if isstruct(filtered), return, end

filtered = filter_by_num(qset, 'Patient.SurvivalTime', ...
	'survival|survival\s*time', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_num(qset, 'Patient.Censored', ...
	'censored|censoring', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_num(qset, 'Patient.AgeAtLastFollowup', ...
	'age\s*at\s*last\s*followup|age\s*at\s*followup', predicate);
if isstruct(filtered), return, end

filtered = filter_by_num(qset, 'Patient.AgeAtDeath', ...
	'age\s*at\s*death|death\s*age', predicate);
if isstruct(filtered), return, end

filtered = filter_by_num(qset, 'Patient.AgeAtDiagnosis', ...
	'age\s*at\s*diagnosis|diagnosis\s*age|age\s*at\s*dx|dx\s*age', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_date(qset, 'Patient.DateOfBirth', ...
	'date\s*of\s*birth|birth\s*date|birth|born', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_date(qset, 'Patient.DateOfDeath', ...
	'date\s*of\s*death|death\s*date|death|died', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_date(qset, 'Patient.DateOfDiagnosis', ...
	'date\s*of\s*diagnosis|diagnosis\s*date|diagnosed', predicate);
if isstruct(filtered), return, end

filtered = filter_by_date(qset, 'Patient.DateOfLastFollowup', ...
	'date\s*of\s*(last)?\s*followup|last\s*followup', predicate);
if isstruct(filtered), return, end
	






filtered = filter_by_string(qset, 'Sample.ID', 'sample|sample\s*id', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_string(qset, 'Sample.ReplicateID', ...
	'replicate\s*id|replicate', predicate);
if isstruct(filtered), return, end

filtered = filter_by_string(qset, 'Sample.Type', 'sample\s*type', predicate);
if isstruct(filtered), return, end

filtered = filter_by_string(qset, 'Sample.Filename', ...
	'filename|file', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_string(qset, 'Sample.Channel', ...
	'channel|sample\s*channel', predicate);
if isstruct(filtered), return, end
	
filtered = filter_by_string(qset, 'Sample.Organ', ...
	'organ|body\s*part|tissue', predicate);
if isstruct(filtered), return, end

	

error 'Unrecognized predicate.';
filtered = [];
return;








function filtered = filter_by_string(qset, field, field_aliases, predicate)

filtered = [];
tokens = regexpi(predicate, ...
	['^\s*(' field_aliases ')\s*=+\s*(.*?)\s*$'], 'tokens');
if length(tokens) == 1
	tokens = tokens{1}; cmp = tokens{2};
	cmp = regexprep(cmp, '(N/A|unknown)', '-', 'ignorecase');
	eval(['filtered = filter_struct(qset, strcmpi(cmp, qset.' field '));']);
	return;
end

tokens = regexpi(predicate, ...
	['^\s*(' field_aliases ')\s*[~!]=\s*(.*?)\s*$'], 'tokens');
if length(tokens) == 1
	tokens = tokens{1}; cmp = tokens{2};
	cmp = regexprep(cmp, '(N/A|unknown)', '-', 'ignorecase');
	eval(['filtered = filter_struct(qset, ~strcmpi(cmp, qset.' field '));']);
	return;
end

tokens = regexpi(predicate, ...
	['^\s*(' field_aliases ')\s*~+\s*(.*?)\s*$'], 'tokens');
if length(tokens) == 1
	tokens = tokens{1}; regex = tokens{2};
	eval(['f = qset.' field ';']);
	selected = false(length(f), 1);
	for k = 1:length(f)
		if regexpi(f{k}, regex), selected(k) = true; end
	end
	filtered = filter_struct(qset, selected);
	return;
end

tokens = regexpi(predicate, ...
	['^\s*(' field_aliases ')\s*!~+\s*(.*?)\s*$'], 'tokens');
if length(tokens) == 1
	tokens = tokens{1}; regex = tokens{2};
	eval(['f = qset.' field ';']);
	selected = false(length(f), 1);
	for k = 1:length(f)
		if regexpi(f{k}, regex), selected(k) = true; end
	end
	filtered = filter_struct(qset, ~selected);
	return;
end

return;



function filtered = filter_by_num(qset, field, field_aliases, predicate)

filtered = [];
tokens = regexp(predicate, ...
	['^\s*(' field_aliases ')\s*=+\s*(.+?)\s*$'], 'tokens');
if length(tokens) == 1
	tokens = tokens{1}; cmp = tokens{2};
	eval(['filtered = filter_struct(qset, qset.' field ' == ' cmp ');']);
	return;
end

tokens = regexp(predicate, ...
	['^\s*(' field_aliases ')\s*[~!]=+\s*(.+?)\s*$'], 'tokens');
if length(tokens) == 1
	tokens = tokens{1}; cmp = tokens{3};
	eval(['filtered = filter_struct(qset, qset.' field ' ~= ' cmp ');']);
	return;
end

tokens = regexp(predicate, ...
	['^\s*(' field_aliases ')\s*(<|<=|>|>=)\s*(.+?)\s*$'], 'tokens');
if length(tokens) == 1
	tokens = tokens{1}; cmp = tokens{3};
	eval(['filtered = filter_struct(qset, qset.' field ' ' tokens{2} ' ' ...
	      cmp ');']);
	return;
end

return;



function filtered = filter_by_date(qset, field, field_aliases, predicate)

filtered = [];
tokens = regexp(predicate, ...
	['^\s*(' field_aliases ')\s*=+\s*(.+?)\s*$'], 'tokens');
if length(tokens) == 1
	tokens = tokens{1}; cmp = tokens{2};
	
	if strcmp(cmp, '-') || strcmpi(cmp, 'nan')
		eval(['filtered = filter_struct(qset, isnan(qset.' field '));']);
		return;
	end
	
	try cmp = num2str(datenum(cmp));
	catch exception
		error 'Right hand argument uses an unrecognized date format.';
	end
	
	eval(['filtered = filter_struct(qset, qset.' field ' == ' cmp ');']);
	return;
end

tokens = regexp(predicate, ...
	['^\s*(' field_aliases ')\s*(~=|!=)\s*(.+?)\s*$'], 'tokens');
if length(tokens) == 1
	tokens = tokens{1}; cmp = tokens{3};
	
	if strcmp(cmp, '-') || strcmpi(cmp, 'nan')
		eval(['filtered = filter_struct(qset, ~isnan(qset.' field '));']);
		return;
	end
	
	try cmp = num2str(datenum(cmp));
	catch exception
		error 'Right hand argument uses an unrecognized date format.';
	end

	eval(['filtered = filter_struct(qset, qset.' field ' ~= ' cmp ');']);
	return;
end

tokens = regexp(predicate, ...
	['^\s*(' field_aliases ')\s*(<|<=|>|>=)\s*(.+?)\s*$'], 'tokens');
if length(tokens) == 1
	tokens = tokens{1}; cmp = tokens{3};
	
	try cmp = num2str(datenum(cmp));
	catch exception
		error 'Right hand argument uses an unrecognized date format.';
	end
	
	eval(['filtered = filter_struct(qset, qset.' field ' ' tokens{2} ' ' ...
		cmp ');']);
	return;
end

return;


