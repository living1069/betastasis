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
	
	tokens = regexpi(line, '!Sample_characteristics_ch(\d)\s+"(.+?):', ...
		'tokens');
	if length(tokens) == 1
		token = tokens{1};
		channel = token{1};
		ch = token{2};
		
		if regexpi(ch, 'sample.*id'), field = 'SampleID';
		elseif regexpi(ch, 'tissue'), field = 'Tissue';
		elseif regexpi(ch, 'disease.*status'), field = 'DiseaseStatus';
		elseif regexpi(ch, 'tumor.*type'), field = 'TumorType';
		elseif regexpi(ch, 'pathological.*stage'), field = 'PathologicalStage';
		elseif regexpi(ch, 'sample.*source'), field = 'SampleSource';
		else continue; end
			
		tokens = regexpi(line, '"(.*?)"', 'tokens');
		eval(['meta.' field 'Ch' channel ' = cell(length(tokens), 1);']);
		for k = 1:length(tokens)
			token = tokens{k}; token = token{1};
			t = regexpi(token, [ch ':\s*(.+)'], 'tokens');
			if length(t) == 1
				t = t{1}; t = t{1};
				eval(['meta.' field 'Ch' channel '{k} = t;']);
			else
				eval(['meta.' field 'Ch' channel '{k} = ''-'';']);
			end
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
	
	if regexpi(line, '!Sample_description\s+');
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
			
			% Convert the description field name to a valid Matlab field name.
			field = strrep(key, ' ', '_');
			field = strrep(field, '?', '');
			field = strrep(field, '/', '_');
			field = strrep(field, '=', '');
			field = strrep(field, '#', 'num');
			field = strrep(field, '(', '');
			field = strrep(field, ')', '');
			
			if length(field) > namelengthmax, continue, end
			
			if ~isfield(meta, field)
				eval(['meta.' field ' = repmat({''-''}, G, 1);']);
			end
			eval(['meta.' field '{k} = val;']);
		end
		
		continue;
	end
	
	if regexpi(line, '!series_matrix_table_begin')
		break;
	end

end

fclose(fid);

