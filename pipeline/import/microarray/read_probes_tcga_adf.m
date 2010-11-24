
% This function extracts microarray probe information from TCGA microarray
% definition files.
% 
% Inputs:
%     adf_file - Filesystem path to the ADF file that contains probe metadata.
%     fasta_file - Filesystem path to the FASTA file that contains the probe
%         sequences. This argument can be left out if the ADF file contains
%         the probe sequences.
%
% Outputs:
%     probes - Data structure that contains the probe positions and sequences.
%
% Author: Matti Annala <matti.annala@tut.fi>

function probes = read_probes_tcga_adf(adf_file, fasta_file)

if nargin >= 2
	probe_seqs = cell(10, 1);
	probe_id_to_idx = containers.Map();
					
	fprintf(1, 'Probe sequences read: ');

	probes_loaded = 0;
	while 1
		[headers, sequences] = fastaread(fasta_file, 'blockread', ...
			[probes_loaded + 1, probes_loaded + 1000]);
		
		probe_seqs((probes_loaded+1):(probes_loaded+length(sequences))) = ...
			sequences;
			
		for k = 1:length(headers)
			probe_id_to_idx(headers{k}) = probes_loaded + k;
		end
		
		for k = 1:floor(log10(probes_loaded))+1
			fprintf(1, '\b');
		end

		probes_loaded = probes_loaded + length(sequences);
		fprintf(1, '%d', probes_loaded);
		
		if length(sequences) < 1000, break, end
	end

	fprintf(1, '\n');
end



fid = fopen(adf_file);
header = fgetl(fid);
fseek(fid, 0, -1);

posx_column = 0;
posy_column = 0;
probe_id_column = 0;
probe_seq_column = 0;

% Figure out which columns we should read.
column_names = regexp(header, '(\w+)(\t|$)', 'tokens');
column_order = [];
for k = 1:length(column_names)
	column_name = column_names{k};
	column_name = column_name{1};
	
	if regexpi(column_name, '^x$')
		posx_column = k;
		column_order(end + 1) = 1;
	elseif regexpi(column_name, '^col$')
		posx_column = k;
		column_order(end + 1) = 1;
	elseif regexpi(column_name, '^y$')
		posy_column = k;
		column_order(end + 1) = 2;
	elseif regexpi(column_name, '^row$')
		posy_column = k;
		column_order(end + 1) = 2;
	elseif regexpi(column_name, '^reporter(_id)?$')
		probe_id_column = k;
		column_order(end + 1) = 3;
	elseif regexpi(column_name, '^ReportSequence$')
		probe_seq_column = k;
		column_order(end + 1) = 4;
	elseif regexpi(column_name, '^ProbeName$')
		probe_id_column = k;
		column_order(end + 1) = 3;
	end
end

if posx_column == 0 || posy_column == 0 || probe_id_column == 0
	error 'Invalid file format, vital columns missing.';
end

if probe_seq_column == 0 && nargin < 2
	error(['ADF file does not contain probe sequences and a FASTA file was ' ...
	       'not specified as an argument.']);
end

fprintf(1, 'Detected file format with %d columns.\n', length(column_names));
fprintf(1, 'Probe X and Y positions stored in columns %d and %d.\n', ...
	posx_column, posy_column);
fprintf(1, 'Probe IDs stored in column %d.\n', probe_id_column);

if probe_seq_column ~= 0
	fprintf(1, 'Probe sequences stored in column %d.\n', probe_seq_column);
end

% Construct a format string for parsing.
parse_format = '';
for k = 1:length(column_names)
	if k == posx_column || k == posy_column
		parse_format = [parse_format '%d'];
	elseif k == probe_id_column || k == probe_seq_column
		parse_format = [parse_format '%s'];
	else
		parse_format = [parse_format '%*s'];
	end
end

data = textscan(fid, parse_format, 'HeaderLines', 1, 'Delimiter', '\t', ...
	'BufSize', 16384);
fclose(fid);

fprintf(1, 'Found %d probes in ADF file.\n', length(data{1}));

probe_x = data{find(column_order == 1)};
probe_y = data{find(column_order == 2)};
probe_ids = data{find(column_order == 3)};

if nargin < 2
	probe_seq = upper(data{find(column_order == 4)});
	valid = ~strcmp('', probe_seq);
	probes.Sequence = probe_seq(valid);
else
	valid = find(probe_id_to_idx.isKey(probe_ids));
	fprintf(1, '%d (%d%%) of the probes had an associated sequence.\n', ...
		length(valid), round(length(valid) / length(probe_ids) * 100));

	indices = cell2mat(probe_id_to_idx.values(probe_ids(valid)));
	
	probes.Sequence = upper(probe_seqs(indices));
end

%probes.ID = probe_ids(valid);
probes.XPos = probe_x(valid);
probes.YPos = probe_y(valid);

