
% READ_NCBI_REFSEQ   Parse an NCBI RefSeq transcriptome build into memory
% 
%    REFSEQ = READ_NCBI_REFSEQ(FILEPATH) parses the compiled RefSeq
%    transcriptome build stored at FILEPATH. Genes, transcripts and exons are
%    read from the transcriptome build, and are organized into the Matlab
%    data structure REFSEQ.
%
%    Only transcriptomic information is read from the file: genomic coordinates
%    and other DNA level information must be read from another source, for
%    instance using READ_UCSC_REFGENE.

% Author: Matti Annala <matti.annala@tut.fi>

function refseq = read_ncbi_refseq(filepath)

genes = struct;
genes.Name = cell(50000, 1);
genes.FullName = cell(50000, 1);
genes.Synonyms = cell(50000, 1);
genes.EntrezID = cell(50000, 1);
genes.TranscriptCount = zeros(50000, 1);
genes.Transcripts = zeros(50000, 10);

transcripts = struct;
transcripts.Name = cell(100000, 1);
transcripts.Sequence = cell(100000, 1);
transcripts.Gene = nan(100000, 1);
transcripts.CDS = nan(100000, 2);
transcripts.Exons = cell(100000, 1);
transcripts.ExonPos = cell(100000, 1);

exons = struct;
exons.ID = cell(500000, 1);

exon_pos = nan(500000, 2);


gene_count = 0;
transcript_count = 0;
exon_count = 0;

gene_map = containers.Map;
transcript_map = containers.Map;


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
			
			tokens = regexp(line, '/note="(.+)"', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1};
				genes.FullName{gene_idx} = tokens{1};
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

fprintf(1, '\n');
fclose(rna_gbff);

% Squeeze any preallocated space out of the data structures.
genes.Name = genes.Name(1:gene_count);
genes.FullName = genes.FullName(1:gene_count);
genes.Synonyms = genes.Synonyms(1:gene_count);
genes.EntrezID = genes.EntrezID(1:gene_count);
genes.TranscriptCount = genes.TranscriptCount(1:gene_count);
genes.Transcripts = genes.Transcripts(1:gene_count, :);

transcripts.Name = transcripts.Name(1:transcript_count);
transcripts.Sequence = transcripts.Sequence(1:transcript_count);
transcripts.Gene = transcripts.Gene(1:transcript_count);
transcripts.CDS = transcripts.CDS(1:transcript_count, :);
transcripts.Exons = transcripts.Exons(1:transcript_count);
transcripts.ExonPos = transcripts.ExonPos(1:transcript_count);

exons.ID = exons.ID(1:exon_count);
exon_pos = exon_pos(1:exon_count, :);

% Assign exon sequences.
exons.Sequence = cell(length(exons.ID), 1);
for t = 1:transcript_count
	ex = transcripts.Exons{t};
	for e = 1:length(ex)
		idx = ex(e);
		exons.Sequence{idx} = transcripts.Sequence{t}( ...
			exon_pos(idx, 1):exon_pos(idx, 2));
	end
end

% Construct associations from genes to transcripts.
for k = 1:length(transcripts.Gene)
	idx = transcripts.Gene(k);
	genes.Transcripts(idx, genes.TranscriptCount(idx) + 1) = k;
	genes.TranscriptCount(idx) = genes.TranscriptCount(idx) + 1;
end

% Sort the genes in alphabetical order.
[genes.Name, order] = sort(genes.Name);
genes.FullName = genes.FullName(order);
genes.Synonyms = genes.Synonyms(order);
genes.EntrezID = genes.EntrezID(order);
genes.TranscriptCount = genes.TranscriptCount(order);
genes.Transcripts = genes.Transcripts(order, :);

% We must also remember to fix the gene references in transcript structures.
inv_order = 1:length(order);
inv_order(order) = inv_order;

transcripts.Gene = inv_order(transcripts.Gene)';

% Convert Entrez IDs into integers.
entrez = nan(length(genes.EntrezID), 1);
for g = 1:length(genes.EntrezID)
	entrez(g) = str2double(genes.EntrezID{g});
end
genes.EntrezID = entrez;

% Link exons with genes.
exons.Gene = nan(length(exons.ID), 1);
for t = 1:transcript_count
	tx_exons = transcripts.Exons{t};
	exons.Gene(tx_exons) = transcripts.Gene(t);
end

% Remove duplicate exons with identical identifiers. We also perform a sanity
% check where we make sure that exons with an identical identifier really
% have identical sequences.
keep_exons = [];
exon_map = nan(length(exons.Gene), 1);
for g = 1:length(genes.Name)
	gene_exons = find(exons.Gene == g);
	if length(gene_exons) == 0, continue, end
	
	gex_ids = exons.ID(gene_exons);
	uniq_gex_ids = unique(gex_ids);
	
	for e = 1:length(uniq_gex_ids)
		u = gene_exons(strcmp(uniq_gex_ids{e}, gex_ids));
		seqs = exons.Sequence(u);
		uniq_seqs = unique(seqs);
		if length(uniq_seqs) ~= 1
			% Exons with the same name did not have identical sequences.
			% We add a distinguishing suffix to the exon names.
			for q = 1:length(uniq_seqs)
				seq_class = u(strcmp(uniq_seqs{q}, seqs));
				keep_exons(end+1) = seq_class(1);
				exon_map(seq_class) = length(keep_exons);
				exons.ID{seq_class(1)} = sprintf('%s(%d)', uniq_gex_ids{e}, q);
			end
		else
			keep_exons(end+1) = u(1);
			exon_map(u) = length(keep_exons);
		end
	end
end

exons.ID = exons.ID(keep_exons);
exons.Gene = exons.Gene(keep_exons);
exons.Sequence = exons.Sequence(keep_exons, :);

% Since we dropped out some exons, we must remember to update the exon
% references in the transcript data structures.
for t = 1:length(transcripts.Name)
	transcripts.Exons{t} = exon_map(transcripts.Exons{t});
end

% Finally, we sort the exons based on their name.

	
% Compile the genes, transcripts and exons into one data structure.
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

