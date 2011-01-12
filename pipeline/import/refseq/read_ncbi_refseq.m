
% READ_NCBI_REFSEQ   Parse an NCBI RefSeq transcriptome build into memory
% 
%    REFSEQ = READ_NCBI_REFSEQ(FILEPATH) parses the compiled RefSeq
%    transcriptome build stored at FILEPATH. Genes, transcripts and exons are
%    read from the transcriptome build, and are organized into the Matlab
%    data structure REFSEQ.

% Author: Matti Annala <matti.annala@tut.fi>

function refseq = read_ncbi_refseq(filepath)

genes = struct;
genes.Name = cell(50000, 1);
genes.Synonyms = cell(50000, 1);
genes.EntrezID = cell(50000, 1);
genes.TranscriptCount = NaN(50000, 1);
genes.Transcripts = NaN(50000, 10);

transcripts = struct;
transcripts.Name = cell(100000, 1);
transcripts.Sequence = cell(100000, 1);
transcripts.Gene = NaN(100000, 1);
transcripts.CDS = NaN(100000, 2);

exons = struct;
exons.ID = cell(500000, 1);
exons.Gene = nan(500000, 1);
exons.Position = nan(500000, 2);

gene_count = 0;
transcript_count = 0;
exon_count = 0;

gene_map = containers.Map;
transcript_map = containers.Map;
exon_map = containers.Map;


rna_gbff = fopen(filepath);
if rna_gbff == -1
	fprintf(1, 'File %s could not be opened.', filepath);
	return;
end

fprintf(1, 'Reading transcript annotations into memory...\n');
fprintf(1, 'Transcripts read: 0');

progress_len = 1;
parse_mode = 0;

while 1
	line = fgetl(rna_gbff);
	if line == -1, break, end
	
	while 1
		if parse_mode == 0
			tokens = regexp(line, '^ACCESSION\s+(\w+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; name = tokens{1};
				if transcript_map.isKey(name)
					error('Transcript %s found twice.', name);
				end
				
				transcript_count = transcript_count + 1;
				transcript_map(name) = transcript_count;
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
				tokens = tokens{1};
				exon_start = str2double(tokens{1});
				exon_end = str2double(tokens{2});
				
				exon_count = exon_count + 1;
				exons.Gene(exon_count) = transcripts.Gene(transcript_count);
				exons.Position(exon_count, :) = [exon_start exon_end];
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
			tokens = regexp(line, '/number="(.+)"', 'tokens');
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

fprintf(1, '\n');
fclose(rna_gbff);

genes.Name = genes.Name(1:gene_count);
genes.EntrezID = genes.EntrezID(1:gene_count);
genes.TranscriptCount = genes.TranscriptCount(1:gene_count);
genes.Transcripts = genes.Transcripts(1:gene_count, :);

transcripts.Name = transcripts.Name(1:transcript_count);
transcripts.Sequence = transcripts.Sequence(1:transcript_count);
transcripts.Gene = transcripts.Gene(1:transcript_count);
transcripts.CDS = transcripts.CDS(1:transcript_count, :);

exons.ID = exons.ID(1:exon_count);
exons.Gene = exons.Gene(1:exon_count);
exons.Position = exons.Position(1:exon_count, :);

for k = 1:length(transcripts.Gene)
	idx = transcripts.Gene(k);
	genes.Transcripts(idx, genes.TranscriptCount(idx) + 1) = k;
	genes.TranscriptCount(idx) = genes.TranscriptCount(idx) + 1;
end

% Sort the genes in alphabetical order.
[genes.Name, order] = sort(genes.Name);
genes.EntrezID = genes.EntrezID(order);
genes.TranscriptCount = genes.TranscriptCount(order);
genes.Transcripts = genes.Transcripts(order, :);

% We must also remember to fix the gene references in exon and transcript
% structures.
inv_order = 1:length(order);
inv_order(order) = inv_order;

transcripts.Gene = inv_order(transcripts.Gene);
exons.Gene = inv_order(exons.Gene);

refseq.Genes = genes;
refseq.Transcripts = transcripts;
refseq.Exons = exons;

	




	


function sequence = read_transcript_sequence(rna_gbff)

sequence = '';

while 1
	line = fgetl(rna_gbff);
	if line == -1, break, end
	
	% The transcript sequence is composed of the normal nucleotides TCGA and
	% the ambiguous nucleotide codes specified by IUPAC.
	matches = regexp(line, '[tcgarykmswbdhvn]+', 'match');
	if length(matches) > 0
		sequence = strcat(sequence, matches{:});
	else
		break;
	end
end

