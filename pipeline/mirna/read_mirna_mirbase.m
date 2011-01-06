function [pre_mirnas, mirnas] = read_mirna_mirbase(filepath, organism)

mirna_dat = fopen(filepath);
if mirna_dat == -1
	fprintf(1, 'File %s could not be opened.', filepath);
	return;
end

pre_mirnas = struct;
pre_mirnas.Accession = cell(10, 1);
pre_mirnas.Name = cell(10, 1);
pre_mirnas.Sequence = cell(10, 1);
pre_mirnas.MatureCount = zeros(10, 1);
pre_mirnas.Matures = zeros(10, 1);
pre_mirnas.MatureOffsets = zeros(10, 1);
pre_mirnas.EntrezGene = zeros(10, 1);

mature_ends = zeros(10, 1);

mirnas = struct;
mirnas.Name = cell(10, 1);
mirnas.Accession = cell(10, 1);
mirnas.Sequence = cell(10, 1);

p = 1;
m = 1;

fprintf(1, 'Reading miRNA annotations into memory...\n');
fprintf(1, 'pre-miRNA entries read: 0');

progress_len = 1;

mirna_to_idx = containers.Map;

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
				pre_mirnas.Name{p} = tokens{1};
			end
		else
			pre_mirnas.Name{p} = tokens{1};
		end
		continue;
	end
		
	matches = regexp(line, '^AC   (\S+);', 'tokens');
	if length(matches) == 1
		tokens = matches{1};
		pre_mirnas.Accession{p} = tokens{1};
		pre_mirnas.EntrezGene(p) = NaN;
		continue;
	end
	
	matches = regexp(line, '^DR   ENTREZGENE; (\d+)', 'tokens');
	if length(matches) == 1
		tokens = matches{1};
		pre_mirnas.EntrezGene(p) = str2double(tokens{1});
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
					pre_mirnas.Name{p});
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
				mirnas.Name{m} = name;
				mirnas.Sequence{m} = '';
				mirnas.Accession{m} = mirna_accession;
				m = m + 1;
			end
			
			pre_mirnas.MatureCount(p) = pre_mirnas.MatureCount(p) + 1;
			pre_mirnas.Matures(p, pre_mirnas.MatureCount(p)) = idx;

			tokens = matches{1};
			pre_mirnas.MatureOffsets(p, pre_mirnas.MatureCount(p)) = ...
				str2double(tokens{1});
			mature_ends(p, pre_mirnas.MatureCount(p)) = str2double(tokens{2});
			break;
		end
	end
	
	if regexp(line, '^SQ   Sequence')
		pre_mirnas.Sequence{p} = read_sequence(mirna_dat);
		preseq = pre_mirnas.Sequence{p};
		
		for k = 1:pre_mirnas.MatureCount(p)
			idx = pre_mirnas.Matures(p, k);
			if mature_ends(p, k) > length(preseq)
				fprintf(1, 'Mature miRNA %s extends beyond pre-miRNA:\n', ...
					mirnas.Name{idx});
			end

			seq = preseq(pre_mirnas.MatureOffsets(p, k):mature_ends(p, k));
			if length(mirnas.Sequence{idx}) > 0 && ...
				~strcmp(seq, mirnas.Sequence{idx})
				error('Mismatching mature sequences for miRNA %s.', ...
					mirnas.Name{idx});
			end
			
			mirnas.Sequence{idx} = seq;
		end
		
		p = p + 1;
		pre_mirnas.MatureCount(p) = 0;
		
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
