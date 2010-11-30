function meta = meta_geo(qset, geo_series_matrix)

meta = qset;
geo = read_geo_series_matrix(geo_series_matrix)

if isfield(meta.Sample, 'Filename')
	S = length(meta.Sample.Filename);
end

geo_to_row = containers.Map(geo.SampleAccession, ...
	num2cell(1:length(geo.SampleAccession)));;
sample_to_row = zeros(S, 1);

if isfield(meta.Sample, 'Filename') && regexp(meta.Sample.Filename{1}, 'GSM\d+')
	for s = 1:length(meta.Sample.Filename)
		tokens = regexp(meta.Sample.Filename{s}, '(GSM\d+)', 'tokens');
		if length(tokens) ~= 1, continue, end
		
		token = tokens{1}; accession = token{1};
		sample_to_row(s) = geo_to_row(accession);
	end
else
	error 'Could not find a link between samples and the GEO series matrix.';
end

found = find(sample_to_row ~= 0);
found_rows = sample_to_row(found);

if isfield(geo, 'SampleIDCh1')
	meta.Sample.ID(found) = geo.SampleIDCh1(found_rows);
end

if isfield(geo, 'TumorTypeCh1')
	meta.Sample.Type(found) = geo.TumorTypeCh1(found_rows);
elseif isfield(geo, 'DiseaseStatusCh1')
	meta.Sample.Type(found) = geo.DiseaseStatusCh1(found_rows);
elseif isfield(geo, 'TissueCh1')
	meta.Sample.Type(found) = geo.TissueCh1(found_rows);
end

