
% READ_GEO_SERIES_MATRIX      Reads data from a GEO series matrix
%
%   [GEO, INFO] = READ_GEO_SERIES_MATRIX(PATH)

% Author: Matti Annala <matti.annala@tut.fi>

function [geo, info] = read_geo_series_matrix(filepath, varargin)

read_data = false;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'data')
		read_data = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

geo = struct;
meta = struct;
info = struct;

fid = fopen(filepath);

while 1
	line = fgetl(fid);
	if line == -1, break; end
	
	tokens = regexpi(line, '!Series_title\s+"(.+)"', 'tokens');
	if length(tokens) == 1
		tokens = tokens{1}; info.title = tokens{1};
		continue;
	end
	
	tokens = regexpi(line, '!Series_geo_accession\s+"(.+)"', 'tokens');
	if length(tokens) == 1
		tokens = tokens{1}; info.series_accession = tokens{1};
		continue;
	end

	tokens = regexpi(line, '!Series_submission_date\s+"(.+)"', 'tokens');
	if length(tokens) == 1
		tokens = tokens{1}; info.submission_date = tokens{1};
		continue;
	end

	tokens = regexpi(line, '!Series_last_update_date\s+"(.+)"', 'tokens');
	if length(tokens) == 1
		tokens = tokens{1}; info.last_update_date = tokens{1};
		continue;
	end

	if regexpi(line, '!Sample_title\s+');
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		meta.sample_title = cell(1, length(tokens));
		for k = 1:length(tokens)
			token = tokens{k}; meta.sample_title{k} = token{1};
		end
		continue;
	end

	if regexpi(line, '!Sample_geo_accession\s+');
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		meta.sample_accession = cell(1, length(tokens));
		for k = 1:length(tokens)
			token = tokens{k}; meta.sample_accession{k} = token{1};
		end
		continue;
	end
	
	tokens = regexpi(line, '!Sample_organism_ch(\d)', 'tokens');
	if length(tokens) == 1
		token = tokens{1};
		org = token{1};
		
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		eval(['meta.organism_ch' channel ' = cell(1, length(tokens));']);
		for k = 1:length(tokens)
			token = tokens{k}; eval(['meta.organism_ch' channel ...
				'{k} = token{1};']);
		end
		continue;
	end
	
	tokens = regexpi(line, '!Sample_characteristics_ch(\d)\s+"[^"]*:','tokens');
	if length(tokens) == 1
		token = tokens{1}; channel = token{1};
		G = length(meta.sample_accession);
		
		tokens = regexpi(line, '"(.*?)"', 'tokens');
		if length(tokens) ~= G
			error 'Invalid number of sample description columns.';
		end
		
		for k = 1:length(tokens)
			token = tokens{k}; token = token{1};
			t = regexpi(token, '(.+?):\s*(.+)', 'tokens');
			if length(t) ~= 1, continue, end
			
			t = t{1}; key = t{1}; val = t{2};
			
			field = valid_field_name(key);
			field = [field '_ch' channel];
			
			if ~isfield(meta, field)
				eval(['meta.' field ' = repmat({''''}, 1, G);']);
			end
			eval(['meta.' field '{k} = val;']);
		end
		continue;
	end

	tokens = regexpi(line, '!Sample_label_ch(\d)', 'tokens');
	if length(tokens) == 1
		token = tokens{1}; channel = token{1};
		
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		eval(['meta.label_ch' channel ' = cell(1, length(tokens));']);
		for k = 1:length(tokens)
			token = tokens{k}; eval(['meta.label_ch' channel ...
				'{k} = token{1};']);
		end
		continue;
	end
	
	tokens = regexpi(line, '!Sample_source_name_ch(\d)', 'tokens');
	if length(tokens) == 1
		token = tokens{1}; channel = token{1};
		
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		eval(['meta.source_name_ch' channel ' = cell(1, length(tokens));']);
		for k = 1:length(tokens)
			token = tokens{k}; eval(['meta.source_name_ch' channel ...
				'{k} = token{1};']);
		end
		continue;
	end
	
	if regexpi(line, '!Sample_description\s+"[^"]*:');
		G = length(meta.sample_accession);
		
		tokens = regexpi(line, '"(.*?)"', 'tokens');
		if length(tokens) ~= G
			error 'Invalid number of sample description columns.';
		end
		
		for k = 1:length(tokens)
			token = tokens{k}; token = token{1};
			t = regexpi(token, '(.+?):\s*(.+)', 'tokens');
			if length(t) ~= 1, continue, end
				
			
			t = t{1}; key = t{1}; val = t{2};
			
			field = valid_field_name(key);
			
			if ~isfield(meta, field)
				eval(['meta.' field ' = repmat({''''}, 1, G);']);
			end
			eval(['meta.' field '{k} = val;']);
		end
		continue;
	end
	
	if regexpi(line, '!Sample_description\s+');
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		meta.sample_description = cell(1, length(tokens));
		for k = 1:length(tokens)
			token = tokens{k}; meta.sample_description{k} = token{1};
		end
		continue;
	end

	if regexpi(line, '!series_matrix_table_begin')
		break;
	end

end

geo.meta = meta;

if read_data
	[data, headers] = readtable(fid, 'Numeric', '^GSM.*');
	geo.rows.id = strrep(data{1}, '"', '');

	geo.mean = nan(length(data{1}), length(data)-1);
	for s = 1:length(data)-1
		geo.mean(:, s) = str2double(data{s+1});
	end

	geo.rows.id = geo.rows.id(1:end-1);
	geo.mean = geo.mean(1:end-1, :);
end

fclose(fid);




function field = valid_field_name(str)

% Convert the description field name to a valid Matlab field name.
field = lower(str);
field = strrep(field, ' ', '_');
field = strrep(field, '/', '_');
field = strrep(field, '?', '');
field = strrep(field, '=', '');
field = strrep(field, '#', 'num');
field = strrep(field, '(', '');
field = strrep(field, ')', '');
field = strrep(field, ',', '_');

if length(field) > namelengthmax
	field = field(1:namelengthmax-3);
end

