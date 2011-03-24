function meta = read_geo_series_matrix(filepath)

meta = struct;

fid = fopen(filepath);

while 1
	line = fgetl(fid);
	if line == -1, break; end
	
	tokens = regexpi(line, '!Series_title\s+"(.+)"', 'tokens');
	if length(tokens) == 1
		tokens = tokens{1}; meta.Title = tokens{1};
		continue;
	end
	
	tokens = regexpi(line, '!Series_geo_accession\s+"(.+)"', 'tokens');
	if length(tokens) == 1
		tokens = tokens{1}; meta.SeriesAccession = tokens{1};
		continue;
	end

	tokens = regexpi(line, '!Series_submission_date\s+"(.+)"', 'tokens');
	if length(tokens) == 1
		tokens = tokens{1}; meta.SubmissionDate = tokens{1};
		continue;
	end

	tokens = regexpi(line, '!Series_last_update_date\s+"(.+)"', 'tokens');
	if length(tokens) == 1
		tokens = tokens{1}; meta.LastUpdateDate = tokens{1};
		continue;
	end

	if regexpi(line, '!Sample_title\s+');
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		meta.SampleTitle = cell(length(tokens), 1);
		for k = 1:length(tokens)
			token = tokens{k}; meta.SampleTitle{k} = token{1};
		end
		continue;
	end

	if regexpi(line, '!Sample_geo_accession\s+');
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		meta.SampleAccession = cell(length(tokens), 1);
		for k = 1:length(tokens)
			token = tokens{k}; meta.SampleAccession{k} = token{1};
		end
		continue;
	end
	
	tokens = regexpi(line, '!Sample_organism_ch(\d)', 'tokens');
	if length(tokens) == 1
		token = tokens{1};
		org = token{1};
		
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		eval(['meta.OrganismCh' channel ' = cell(length(tokens), 1);']);
		for k = 1:length(tokens)
			token = tokens{k}; eval(['meta.OrganismCh' channel ...
				'{k} = token{1};']);
		end
		continue;
	end
	
	tokens = regexpi(line, '!Sample_characteristics_ch(\d)\s+"[^"]*:','tokens');
	if length(tokens) == 1
		token = tokens{1}; channel = token{1};
		G = length(meta.SampleAccession);
		
		tokens = regexpi(line, '"(.*?)"', 'tokens');
		if length(tokens) ~= G
			error 'Invalid number of sample description columns.';
		end
		
		for k = 1:length(tokens)
			token = tokens{k}; token = token{1};
			t = regexpi(token, '(.+?):\s*(.+)', 'tokens');
			if length(t) ~= 1
				eval(['meta.' field 'Ch' channel '{k} = ''-'';']);
				continue;
			end
			
			t = t{1}; key = t{1}; val = t{2};
			
			field = valid_field_name(key);
			field = [field 'Ch' channel];
			
			if ~isfield(meta, field)
				eval(['meta.' field ' = repmat({''-''}, G, 1);']);
			end
			eval(['meta.' field '{k} = val;']);
		end
		continue;
	end

	tokens = regexpi(line, '!Sample_label_ch(\d)', 'tokens');
	if length(tokens) == 1
		token = tokens{1}; channel = token{1};
		
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		eval(['meta.LabelCh' channel ' = cell(length(tokens), 1);']);
		for k = 1:length(tokens)
			token = tokens{k}; eval(['meta.LabelCh' channel ...
				'{k} = token{1};']);
		end
		continue;
	end
	
	tokens = regexpi(line, '!Sample_source_name_ch(\d)', 'tokens');
	if length(tokens) == 1
		token = tokens{1}; channel = token{1};
		
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		eval(['meta.SourceNameCh' channel ' = cell(length(tokens), 1);']);
		for k = 1:length(tokens)
			token = tokens{k}; eval(['meta.SourceNameCh' channel ...
				'{k} = token{1};']);
		end
		continue;
	end
	
	if regexpi(line, '!Sample_description\s+"[^"]*:');
		G = length(meta.SampleAccession);
		
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
				eval(['meta.' field ' = repmat({''-''}, G, 1);']);
			end
			eval(['meta.' field '{k} = val;']);
		end
		continue;
	end
	
	if regexpi(line, '!Sample_description\s+');
		tokens = regexpi(line, '"(.+?)"', 'tokens');
		meta.SampleDesc = cell(length(tokens), 1);
		for k = 1:length(tokens)
			token = tokens{k}; meta.SampleDesc{k} = token{1};
		end
		continue;
	end

	
	if regexpi(line, '!series_matrix_table_begin')
		break;
	end

end

fclose(fid);




function field = valid_field_name(str)

% Convert the description field name to a valid Matlab field name.
field = strrep(str, ' ', '_');
field = strrep(field, '/', '_');
field = strrep(field, '?', '');
field = strrep(field, '=', '');
field = strrep(field, '#', 'num');
field = strrep(field, '(', '');
field = strrep(field, ')', '');

ws = (field == '_');
field(find(ws)+1) = upper(field(find(ws)+1));
field(1) = upper(field(1));
field = field(~ws);

field = strrep(field, ',', '_');

if length(field) > namelengthmax
	field = field(1:namelengthmax-3);
end

