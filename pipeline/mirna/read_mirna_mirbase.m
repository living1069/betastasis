function [pre_mirnas, mirnas] = read_mirna_mirbase(filepath, organism)

mirna_dat = fopen(filepath);
if mirna_dat == -1
	fprintf(1, 'File %s could not be opened.', filepath);
	return;
end

pre_mirna_accessions = cell(10, 1);
pre_mirna_names = cell(10, 1);
mature_count = zeros(10, 1);
mature_mirnas = zeros(10, 1);
mature_offsets = zeros(10, 1);
mature_ends = zeros(10, 1);
pre_mirna_sequences = cell(10, 1);

mirna_names = cell(10, 1);
mirna_accessions = cell(10, 1);
mirna_sequences = cell(10, 1);

p = 1;
m = 1;

fprintf(1, 'Reading miRNA annotations into memory...\n');
fprintf(1, 'pre-miRNA entries read: 0');

progress_len = 1;

mirna_to_idx = containers.Map();

while 1
	line = fgetl(mirna_dat);
	if line == -1, break, end
	
	matches = regexp(line, '^ID   (\S+)', 'tokens');
	if length(matches) == 1
		tokens = matches{1};
				
		% If the user only wants miRNA's for a specific organism, we have to
		% check whether this miRNA belong to that organism.
		if nargin == 2
			name = tokens{1};
			if ~strcmp(organism, name(1:length(organism)))
				% Wrong organism, so we skip this entry altogether.
				while 1
					skip_line = fgetl(mirna_dat);
					if skip_line == -1, break, end
					if strcmp('//', skip_line(1:2)), break, end
				end
			else
				pre_mirna_names{p} = tokens{1};
			end
		else
			pre_mirna_names{p} = tokens{1};
		end
		continue;
	end
		
	matches = regexp(line, '^AC   (\S+);', 'tokens');
	if length(matches) == 1
		tokens = matches{1};
		pre_mirna_accessions{p} = tokens{1};
		continue;
	end
	
	matches = regexp(line, '^FT   miRNA\s+(\d+)\.\.(\d+)', 'tokens');
	if length(matches) == 1
		mirna_accession = '';
		
		while 1
			prod_line = fgetl(mirna_dat);
			if prod_line == -1, break, end
				
			% Sanity check to prevent situations where we skip over other
			% data while looking for a /product field.
			if length(regexp(prod_line, '^FT   miRNA\s+(\d+)\.\.(\d+)', ...
				'tokens')) == 1
				fprintf(1, 'WARNING: No /product row found for miRNA %s.\n', ...
					pre_mirna_names{p});
				break;
			end

			acc_matches = regexp(prod_line, '^FT.*/accession="(.+)"', 'tokens');
			if length(acc_matches) == 1
				tokens = acc_matches{1};
				mirna_accession = tokens{1};
			end
			
			prod_matches = regexp(prod_line, '^FT.*/product="(.+)"', 'tokens');
			if length(prod_matches) == 0, continue, end
			
			tokens = prod_matches{1};
			name = tokens{1};
			
			% Check if we have already seen this mature miRNA. A mature miRNA
			% can be produced from multiple different pre-miRNA variants, and
			% can thus appear multiple times in the database.
			idx = m;
			if mirna_to_idx.isKey(name)
				idx = mirna_to_idx(name);
			else
				if strcmp('', mirna_accession)
					error 'Found a mature miRNA without an accession ID.';
				end
				
				mirna_to_idx(name) = m;
				mirna_names{m} = name;
				mirna_sequences{m} = '';
				mirna_accessions{m} = mirna_accession;
				m = m + 1;
			end
			
			mature_count(p) = mature_count(p) + 1;
			mature_mirnas(p, mature_count(p)) = idx;

			tokens = matches{1};
			mature_offsets(p, mature_count(p)) = str2double(tokens{1});
			mature_ends(p, mature_count(p)) = str2double(tokens{2});
			break;
		end
	end
	
	if regexp(line, '^SQ   Sequence')
		pre_mirna_sequences{p} = read_sequence(mirna_dat);
		preseq = pre_mirna_sequences{p};
		
		for k = 1:mature_count(p)
			idx = mature_mirnas(p, k);
			if mature_ends(p, k) > length(preseq)
				fprintf(1, 'Mature miRNA %s extends beyond pre-miRNA:\n', ...
					mirna_names{idx});
			end

			seq = preseq(mature_offsets(p, k):mature_ends(p, k));
			if length(mirna_sequences{idx}) > 0 && ...
				~strcmp(seq, mirna_sequences{idx})
				error('Mismatching mature sequences for miRNA %s.', ...
					mirna_names{idx});
			end
			
			mirna_sequences{idx} = seq;
		end
		
		p = p + 1;
		mature_count(p) = 0;
		
		if mod(p, 100) == 0
			for j = 1:progress_len, fprintf(1, '\b'); end
			fprintf(1, '%d', p);
			progress_len = length(int2str(p));
		end
	end
end

for j = 1:progress_len, fprintf(1, '\b'); end
fprintf(1, '%d', p - 1);
fprintf(1, '\n');

fclose(mirna_dat);

pre_mirnas = struct( ...
	'Accession', { pre_mirna_accessions }, ...
	'Name', { pre_mirna_names }, ...
	'Sequence', { pre_mirna_sequences }, ...
	'MatureCount', mature_count, ...
	'Matures', mature_mirnas, ...
	'MatureOffsets', mature_offsets);

mirnas = struct( ...
	'Name', { mirna_names }, ...
	'Accession', { mirna_accessions }, ...
	'Sequence', { mirna_sequences });
	
return;
	
	


function sequence = read_sequence(file)

sequence = '';

while 1
	line = fgetl(file);
	if line == -1, break, end
		
	matches = regexp(line, '[ucgarykmswbdhvn]+', 'match');
	if length(matches) > 0
		sequence = strcat(sequence, matches{:});
	else
		break;
	end
end

sequence = strrep(upper(sequence), 'U', 'T');
return;
