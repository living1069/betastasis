function [genes, transcripts, exons] = read_ncbi_gbs()

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;

chromosomes = struct;
chromosomes.Name = {};
for k = 1:22, chromosomes.Name{k} = num2str(k); end
chromosomes.Name{23} = 'X';
chromosomes.Name{24} = 'Y';
chromosomes.Name{25} = 'M';

transcripts.Name = {};
transcripts.Sequence = {};
transcripts.Gene = [];
transcripts.ExonCount = [];
transcripts.Exons = [];
transcripts.ExonPositions = [];

exons = struct;
exons.Transcript = [];
exons.Coordinates = [];

gene_map = containers.Map;
transcript_map = containers.Map;
exon_map = containers.Map;




gbs = fopen(filepath);
if gbs == -1, error('File %s could not be opened.', filepath); end

fprintf(1, 'Reading transcript annotations into memory...\n');
fprintf(1, 'Transcripts read: 0');

progress_len = 1;
parse_mode = 0;

while 1
	%if transcript_count == 1000 && ~isempty(transcripts.Sequence{1000})
	%	break;
	%end
	
	line = fgetl(rna_gbff);
	if line == -1, break, end
		
	while 1
		if parse_mode == 0
			tokens = regexp(line, '^ACCESSION\s+(\w+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; name = tokens{1};
				%if transcript_map.isKey(name)
				%	error('Transcript %s found twice.', name);
				%end
				
				transcript_count = transcript_count + 1;
				%transcript_map(name) = transcript_count;
				transcripts.Name{transcript_count} = name;
				break;
			end
			
			tokens = regexp(line, '^     gene');
			if length(tokens) == 1
				gene_idx = NaN;
				parse_mode = 1;
				break;
			end
			
			tokens = regexp(line, '^     exon\s+(\d+)\.\.(\d+)', 'tokens');
			if length(tokens) == 1
				token = tokens{1};
				exon_start = str2double(token{1});
				exon_end = str2double(token{2});

				exon_count = exon_count + 1;
				transcripts.Exons{transcript_count}(end+1) = exon_count;
				transcripts.ExonPos{transcript_count}(end+1, :) = ...
					[exon_start exon_end];
				exon_pos(exon_count, :) = [exon_start exon_end];
				
				parse_mode = 2;
				break;
			end
			
			tokens = regexp(line, '\s+CDS\s+(\d+)\.\.(\d+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1};
				cds_start = str2double(tokens{1});
				cds_end = str2double(tokens{2});
				
				if ~isnan(transcripts.CDS(transcript_count, 1))
					error('Transcript %s has two CDS.', ...
						transcripts.Name{transcript_count});
				end
				
				transcripts.CDS(transcript_count, :) = [cds_start cds_end];
				break;
			end

			
			if regexp(line, '^ORIGIN')
				transcripts.Sequence{transcript_count} = ...
					read_transcript_sequence(rna_gbff);
				
				% Update the progress indicator that is shown to the user.
				if mod(transcript_count, 100) == 0
					for j = 1:progress_len, fprintf(1, '\b'); end
					fprintf(1, '%d', transcript_count);
					progress_len = length(int2str(transcript_count));
				end
			end

			break;
			
			
			
			
		% Parsing a gene feature.
		elseif parse_mode == 1
			tokens = regexp(line, '/gene="(.+)"', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; name = tokens{1};
				if ~gene_map.isKey(name)
					gene_count = gene_count + 1;
					gene_map(name) = gene_count;
					genes.Name{gene_count} = name;
				end
				gene_idx = gene_map(name);

				transcripts.Gene(transcript_count) = gene_idx;
				break;
			end
			
			tokens = regexp(line, '/gene_synonym="(.+)"', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1};
				synonyms = textscan(tokens{1}, '%s', 'Delimiter', ';');
				genes.Synonyms{gene_idx} = synonyms{1};
				break;
			end
			
			tokens = regexp(line, '/db_xref="GeneID:(\d+)"', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; geneid = tokens{1};
				genes.EntrezID{gene_idx} = geneid;
				break;
			end
			
			if length(line) < 10 || any(line(1:10) ~= ' ')
				parse_mode = 0;
				continue;
			end
			
			break;
			
			
		
		
		% Parsing an exon feature.
		elseif parse_mode == 2
			tokens = regexp(line, '/number=(\S+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; name = tokens{1};
				exons.ID{exon_count} = name;
				break;
			end
			
			if length(line) < 10 || any(line(1:10) ~= ' ')
				parse_mode = 0; continue;
			end
			
			break;
		end
	end
end

