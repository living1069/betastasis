
% This function imports an NCBI RefSeq transcriptome build into memory as a
% Matlab structure.
% 
% The input file is assumed to contain only transcripts that are relevant to
% the target organism. This means that a freshly downloaded NCBI RefSeq build
% must first be filtered with the script 'compile_human_refseq.py' in order
% to produce a suitable *.rna.gbff file.
%
% Inputs:
%     filepath - Filesystem path to an *.rna.gbff file that contains a full set
%         of transcriptome annotations for a particular organism.
%
% Outputs:
%     transcriptome - Structure that contains transcript-level information
%         about the organism's transcriptome.
%     exons - Structure that contains the positions of all known exons within
%         the organism's transcriptome.
%
% Author: Matti Annala <matti.annala@tut.fi>

function [genes, transcripts, exons] = read_ncbi_refseq(filepath, orgname)

chromosomes = chromosomes_for_organism(orgname);

genes = struct;
genes.Name = cell(50000, 1);
genes.EntrezID = cell(50000, 1);
genes.TranscriptCount = zeros(50000, 1);
genes.Transcripts = zeros(50000, 10);

transcripts = struct;
transcripts.Name = cell(100000, 1);
transcripts.Sequence = cell(100000, 1);
transcripts.Gene = zeros(100000, 1);
transcripts.Chromosome = zeros(100000, 1);
transcripts.CDS = zeros(100000, 2);

exons = struct;
exons.Transcript = zeros(500000, 1);
exons.Position = zeros(500000, 2);

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
			
			tokens = regexp(line, '\s+\/chromosome="(\w+)"', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; chr = tokens{1};
				if strcmpi(chr, 'Unknown') || strcmpi(chr, 'Un')
					chrnum = 0;
				else
					chrnum = find(strcmp(chr, chromosomes));
				end
				
				if length(chrnum) ~= 1, error('Bad chr %s.', chr); end
				transcripts.Chromosome(transcript_count) = chrnum;
				break;
			end

			tokens = regexp(line, '^\s+gene\s+[<>]?\d+\.\.[<>]?\d+');
			if length(tokens) == 1
				gene_idx = 0;
				parse_mode = 1;
				break;
			end
			
			tokens = regexp(line, '\s+exon\s+(\d+)\.\.(\d+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1};
				exon_start = str2double(tokens{1});
				exon_end = str2double(tokens{2});
				
				exon_count = exon_count + 1;
				exons.Transcript(exon_count) = transcript_count;
				exons.Position(exon_count, :) = [exon_start exon_end];
				break;
			end
			
			tokens = regexp(line, '\s+CDS\s+(\d+)\.\.(\d+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1};
				cds_start = str2double(tokens{1});
				cds_end = str2double(tokens{2});
				
				if transcripts.CDS(transcript_count, 1) ~= 0
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
			
			tokens = regexp(line, '/db_xref="GeneID:(\d+)"', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; geneid = tokens{1};
				genes.EntrezID{gene_idx} = geneid;
				break;
			end
			
			if length(line) < 10 || sum(line(1:10) ~= ' ')
				parse_mode = 0;
				continue;
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
transcripts.Chromosome = transcripts.Chromosome(1:transcript_count);
transcripts.CDS = transcripts.CDS(1:transcript_count, :);

% Transcripts that have no CDS must have NaN in their 'CDS' offsets.
transcripts.CDS(transcripts.CDS(:, 1) == 0, :) = NaN;

exons.Transcript = exons.Transcript(1:exon_count);
exons.Position = exons.Position(1:exon_count, :);

for k = 1:length(transcripts.Gene)
	idx = transcripts.Gene(k);
	genes.Transcripts(idx, genes.TranscriptCount(idx) + 1) = k;
	genes.TranscriptCount(idx) = genes.TranscriptCount(idx) + 1;
end

[genes.Name, order] = sort(genes.Name);
genes.EntrezID = genes.EntrezID(order);
genes.TranscriptCount = genes.TranscriptCount(order);
genes.Transcripts = genes.Transcripts(order, :);

return;
	
	


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

return;





function chromosomes = chromosomes_for_organism(orgname)

chromosomes = {};

if strcmpi(orgname, 'homo sapiens')
	for k = 1:22, chromosomes{k} = num2str(k); end
	chromosomes{23} = 'X';
	chromosomes{24} = 'Y';
	chromosomes{25} = 'M';
elseif strcmpi(orgname, 'mus musculus')
	for k = 1:19, chromosomes{k} = num2str(k); end
	chromosomes{20} = 'X';
	chromosomes{21} = 'Y';
	chromosomes{22} = 'M';
else
	error 'Karyotype of specified organism is not known.'; 
end

return;
