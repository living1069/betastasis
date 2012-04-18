function ps = read_probesets_affy_cdf(cdf_file, probes)

cdf = affyread(cdf_file);

ps_ids = { cdf.ProbeSets.Name }';
ps_ids = ps_ids(~strcmp('', ps_ids));

ps = probesetlookup(cdf, ps_ids);

ps_targets = { ps.Identifier }';
ps_targets = regexprep(ps_targets, '

return;


% Affymetrix microarrays have differences in how their probe description files
% are formatted. The right format must first be detected by examining the file's
% first line, which contains the column headers.
parse_format = '';

fid = fopen(filepath);
header = fgetl(fid);
%fseek(fid, 0, -1);

column_names = textscan(header, '%s', -1, 'Delimiter', '\t');
column_names = column_names{1};

probeset_id_col = 0;
probe_x_col = 0;
probe_y_col = 0;
probe_seq_col = 0;

for k = 1:length(column_names)
	if regexpi(column_names{k}, 'probe.?set.?(id|name)')
		probeset_id_col = k;
	elseif regexpi(column_names{k}, 'probe.x(.pos)?')
		probe_x_col = k;
	elseif regexpi(column_names{k}, 'probe.y(.pos)?')
		probe_y_col = k;
	elseif regexpi(column_names{k}, 'probe.sequence')
		probe_seq_col = k;
	end
end

if ~probeset_id_col || ~probe_x_col || ~probe_y_col || ~probe_seq_col
	error 'Invalid file format, vital columns missing.';
end

% Construct a format string for parsing.
parse_format = '';
for k = 1:length(column_names)
	if k == probe_x_col || k == probe_y_col
		parse_format = [parse_format '%d'];
	elseif k == probeset_id_col || k == probe_seq_col
		parse_format = [parse_format '%s'];
	else
		parse_format = [parse_format '%*s'];
	end
end

% Read probe information from the file once the file format has been determined.
data = textscan(fid, parse_format);
fclose(fid);

probes = struct('XPos', data{2}, ...
                'YPos', data{3}, ...
				'Sequence', { lower(data{4}) });

fprintf(1, 'Found %d probes.\n', length(probes.Sequence));

if nargout > 1
	probeset_names = unique(data{1});
	probeset_map = containers.Map(probeset_names, ...
		num2cell(1:length(probeset_names)));
					
	probesets = struct('Name', { probeset_names }, ...
					   'ProbeCount', zeros(length(probeset_names), 1), ...
					   'Probes', zeros(length(probeset_names), 10));

	fprintf(1, 'Reading probesets from file...\n');
	progress = Progress;

	probeset_names = data{1};

	N = length(probeset_names);
	for k = 1:N
		ps = probeset_map(probeset_names{k});
		probenum = probesets.ProbeCount(ps);
		if probenum >= 20, continue, end
		probesets.Probes(ps, probenum + 1) = k;
		probesets.ProbeCount(ps) = probenum + 1;
		
		progress.update(k / N);
	end
end

