function meta = meta_geo(qset, geo_series_matrix)

meta = qset;
geo = read_geo_series_matrix(geo_series_matrix);

S = length(meta.Sample.Filename);

geo_to_row = containers.Map(geo.SampleAccession, ...
	num2cell(1:length(geo.SampleAccession)));;
sample_to_row = zeros(S, 1);

if regexp(meta.Sample.Filename{1}, 'GSM\d+')
	for s = 1:length(meta.Sample.Filename)
		tokens = regexp(meta.Sample.Filename{s}, '(GSM\d+)', 'tokens');
		if length(tokens) ~= 1, continue, end
		
		token = tokens{1}; accession = token{1};
		if geo_to_row.isKey(accession)
			sample_to_row(s) = geo_to_row(accession);
		end
	end
elseif isfield(meta.Sample, 'GeoAccession')
	for s = 1:length(meta.Sample.GeoAccession)
		tokens = regexp(meta.Sample.GeoAccession{s}, '(GSM\d+)', 'tokens');
		if length(tokens) ~= 1, continue, end
		
		token = tokens{1}; accession = token{1};
		if geo_to_row.isKey(accession)
			sample_to_row(s) = geo_to_row(accession);
		end
	end
else
	error 'Could not find a link between samples and the GEO series matrix.';
end

found = (sample_to_row ~= 0);
found_rows = sample_to_row(found);

if ~isfield(meta, 'Misc'), meta.Misc = struct; end	
	
G = length(geo.SampleAccession);

fields = fieldnames(geo);
for k = 1:length(fields)
	f = getfield(geo, fields{k});
	if length(f) ~= G, continue, end
		
	% If a field by this name already exists, we simply add the new values
	% to the field without clearing any existing ones.
	if ~isfield(meta.Misc, fields{k})
		eval(['meta.Misc.' fields{k} ' = repmat({''-''}, S, 1);']);
	end
	eval(['meta.Misc.' fields{k} '(found) = f(found_rows);']);
end

return;




if isfield(geo, 'SampleIDCh1')
	if ~isfield(meta.Sample, 'ID')
		meta.Sample.ID = repmat({'-'}, S, 1);
	end
	meta.Sample.ID(found) = geo.SampleIDCh1(found_rows);
end

if isfield(geo, 'TumorTypeCh1')
	if ~isfield(meta.Sample, 'Type')
		meta.Sample.Type = repmat({'-'}, S, 1);
	end
	meta.Sample.Type(found) = geo.TumorTypeCh1(found_rows);
elseif isfield(geo, 'DiseaseStatusCh1')
	if ~isfield(meta.Sample, 'Type')
		meta.Sample.Type = repmat({'-'}, S, 1);
	end
	meta.Sample.Type(found) = geo.DiseaseStatusCh1(found_rows);
elseif isfield(geo, 'TissueCh1')
	if ~isfield(meta.Sample, 'Type')
		meta.Sample.Type = repmat({'-'}, S, 1);
	end
	meta.Sample.Type(found) = geo.TissueCh1(found_rows);
end

