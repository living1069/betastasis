
% CREATE_GENE_PROBESETS   Builds custom gene expression probesets from probes
% 
%    PROBESETS = CREATE_GENE_PROBESETS(PROBES, ...) constructs customizable
%    gene expression probesets from the microarray probes described by data
%    structure PROBES. The probesets are designed based on the currently
%    selected transcriptome build (see HELP ORGANISM).
%    
%    The probesets are built by aligning the probe sequences against the
%    reference transcriptome. A probe is then assigned to the probeset of
%    gene Z if and only if that probe's sequence aligns to every transcript of
%    gene Z, and if and only if the probe doesn't align to any other gene's
%    transcript.
%
%    CREATE_GENE_PROBESETS(..., 'MaxMismatches', MAX_MM) allows MAX_MM 
%    nucleotide mismatches when aligning the probes against the transcriptome.
%    Default is 0 nucleotide mismatches.
%
%    CREATE_GENE_PROBESETS(..., 'IncludeGenes', GENES) only builds probesets
%    for GENES. Default is to include all genes.
%
%    CREATE_GENE_PROBESETS(..., 'EveryTranscript', false) removes the
%    requirement that a probe is required to align to every transcript of a
%    gene in order to be included in its probeset. The default option is to
%    require this.
%
%    CREATE_GENE_PROBESETS(..., 'MinProbesetSize', NPROBES) discards probesets
%    with less than NPROBES probes. Default is NPROBES = 4.
%
%    CREATE_GENE_PROBESETS(..., 'AllowNonspecific', ALLOW) lets the user specify
%    whether a probe is allowed to align to transcripts from multiple genes.
%    Default is to not allow this (ALLOW = false).
%
%    CREATE_GENE_PROBESETS(..., 'ReportDir', DIR) makes the function write
%    JSON format reports about the probesets generated for each gene.

% Author: Matti Annala <matti.annala@tut.fi>

function probesets = create_gene_probesets(probes, varargin)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;

include_genes = true(length(genes.Name), 1);
require_every_transcript = false;
min_probes_per_probeset = 1;
max_mismatches = 0;
allow_nonspecific_hits = false;
report_dir = '';

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MaxMismatches')
		max_mismatches = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'IncludeGenes')
		idx = gene_idx(varargin{k+1});
		include_genes(:) = false;
		include_genes(idx(~isnan(idx))) = true;
		continue;
	end
	
	if strcmpi(varargin{k}, 'EveryTranscript')
		require_every_transcript = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'MinProbesetSize')
		min_probes_per_probeset = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'AllowNonspecific')
		allow_nonspecific_hits = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'ReportDir')
		report_dir = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

fprintf('Writing probe sequences to a temporary file...\n');
probes_fasta_file = ptemp;
write_seq_fasta(probes, probes_fasta_file);

fprintf('Aligning probes to transcripts using Bowtie...\n');
alignments_file = ptemp;
index = '/data/csb/tools/bowtie-indexes/homo_sapiens/transcripts.refseq_38';
unix( ...
	sprintf('/data/csb/tools/bowtie-0.12.8/bowtie -p4 -f -B1 -v%d -y --all --suppress 5,6,7,8 %s %s > %s', ...
	max_mismatches, index, probes_fasta_file, alignments_file));

fprintf('Constructing gene expression probesets based on alignments...\n');

transcript_map = containers.Map(transcripts.Name, ...
	num2cell(1:length(transcripts.Name)));

probesets = struct;
probesets.gene = genes.Name;
probesets.probecount = zeros(length(genes.Name), 1);
probesets.probes = zeros(length(genes.Name), 1);

% Read the Bowtie alignments into Matlab.
alignments = fopen(alignments_file);
data = textscan(alignments, '%d %*s %s %d');
fclose(alignments);

probe_indices = data{1};
tx_id = data{2};
offsets = data{3};

% Map the transcript names found in the Bowtie alignments file into transcript
% indices. It is faster to map them all in one go, rather than inside the loop.
transcript_indices = cell2mat(transcript_map.values(tx_id));

% Find sequential lines that indicate multiple alignments for one probe.
% These are identified by looking at the probe ID column.
run_ends = [ find(probe_indices(1:end-1) ~= probe_indices(2:end)); ...
	length(probe_indices) ];
run_lengths = diff([0; run_ends]);

if ~isempty(report_dir)
	report = struct;
	report.TxProbePosition = cell(length(transcripts.Name), 1);
	report.TxProbeAccepted = cell(length(transcripts.Name), 1);
	report.TxProbeSeq = cell(length(transcripts.Name), 1);
	
	for k = 1:length(report.TxProbeSeq)
		report.TxProbeSeq{k} = {};
	end
	
	probe_len = zeros(length(probes.Sequence), 1);
	for k = 1:length(probes.Sequence)
		probe_len(k) = length(probes.Sequence{k});
	end
end

progress = Progress;

pos = 1;
for r = 1:length(run_lengths)
	% Determine transcript IDs for the transcripts that the probe matches with.
	probe_target_transcripts = transcript_indices(pos:pos+run_lengths(r)-1);
	
	% Determine a set of genes so that at least one of the target transcripts
	% belongs to each gene in the set.
	probe_target_genes = unique(transcripts.Gene(probe_target_transcripts));
	
	if require_every_transcript
		% For each potential target gene, check if the probe targets every
		% transcript variant of the gene.
		matched_genes = [];
		for g = 1:length(probe_target_genes)
			idx = probe_target_genes(g);
			if isempty(setdiff( ...
				genes.Transcripts(idx, 1:genes.TranscriptCount(idx)), ...
				probe_target_transcripts))
				
				matched_genes(end + 1) = idx;
			end
		end
	else
		matched_genes = probe_target_genes;
	end
	
	probe_accepted = false;
	
	% Don't add the probe to any probeset if it matches with more than one
	% gene. Note that doing this check has a huge positive impact on
	% correlations with RNA-Seq data.
	if allow_nonspecific_hits || length(matched_genes) == 1
		for g = 1:length(matched_genes)
			idx = matched_genes(g);
			probe_num = probesets.probecount(idx);
			if probe_num < 40
				probesets.probes(idx, probe_num + 1) = probe_indices(pos);
				probesets.probecount(idx) = probe_num + 1;
				probe_accepted = true;
			end
		end
	end
	
	% Store the probe target positions for report generation.
	if ~isempty(report_dir)
		for p = pos:pos+run_lengths(r)-1
			tx = transcript_indices(p);
			report.TxProbePosition{tx}(end+1, :) = ...
				[offsets(p), offsets(p) + probe_len(probe_indices(p)) - 1];
			report.TxProbeAccepted{tx}(end+1) = probe_accepted;
			report.TxProbeSeq{tx}{end+1} = probes.Sequence{probe_indices(p)};
		end
	end
	
	progress.update(r / length(run_lengths));
	pos = pos + run_lengths(r);
end

if any(~include_genes)
	probesets.probes(~include_genes, :) = 0;
	probesets.probecount(~include_genes) = 0;
end

if min_probes_per_probeset > 1
	fprintf('Filtering out probesets with less than %d probes...\n', ...
		min_probes_per_probeset);
	probesets.probecount(probesets.probecount < min_probes_per_probeset) = 0;
	probesets.probes(probesets.probecount < min_probes_per_probeset, :) = 0;
end

probesets.type = 'Gene expression';
probesets.organism = organism.Name;

% Remove temporary files.
safe_delete(probes_fasta_file);
safe_delete(alignments_file);

% Finally we write a report about the generated probesets.
if ~isempty(report_dir)
	fprintf(1, 'Writing a probeset report to %s...\n', report_dir);
	
	[~, ~] = mkdir(report_dir);
	
	fid = fopen([report_dir '/genelist.json'], 'w');
	fprintf(fid, '{ "genes": [ ');
	
	for g = 1:length(genes.Name)
		plot_probe_reads(report, g, report_dir);
		if g == 1 
			fprintf(fid, '"%s"', organism.Genes.Name{g});
		else
			fprintf(fid, ', "%s"', organism.Genes.Name{g});
		end
	end
	
	fprintf(fid, ' ] }\n');
	fclose(fid);
end








function [] = plot_probe_reads(report, gene_idx, report_dir)
	
global organism;
transcripts = organism.Transcripts;
exons = organism.Exons;

tx_indices = find(transcripts.Gene == gene_idx);

longest_tx = 0;
for t = 1:length(tx_indices)
	tx_seq = transcripts.Sequence{tx_indices(t)};
	longest_tx = max(longest_tx, length(tx_seq));
end

gene_name = organism.Genes.Name{gene_idx};
[~, ~] = mkdir([report_dir '/' lower(gene_name(1))]);

fid = fopen([report_dir '/' lower(gene_name(1)) '/' ...
	organism.Genes.Name{gene_idx} '.json'], 'W');
fprintf(fid, '{\n');

for t = 1:length(tx_indices)
	tx_idx = tx_indices(t);
	tx_seq = transcripts.Sequence{tx_idx};
	
	probe_pos = report.TxProbePosition{tx_idx};
	accepted = report.TxProbeAccepted{tx_idx};
	probe_seq = report.TxProbeSeq{tx_idx};
	
	fprintf(fid, '"%s": {\n', transcripts.Name{tx_idx});
	fprintf(fid, '\t"length": %d,\n', length(tx_seq));
	fprintf(fid, '\t"exons": [');
	
	tx_exons = transcripts.Exons{tx_idx};
	for m = 1:length(tx_exons)
		exon_id = exons.ID{tx_exons(m)};
		if m == length(tx_exons)
			fprintf(fid, '"%s"],\n', exon_id);
		else
			fprintf(fid, '"%s", ', exon_id);
		end
	end
	
	if length(tx_exons) == 0
		fprintf(fid, '],\n');
	end
	
	fprintf(fid, '\t"exon_pos": [');
	
	tx_exon_pos = transcripts.ExonPos{tx_idx};
	for m = 1:size(tx_exon_pos, 1)
		left = tx_exon_pos(m, 1);
		right = tx_exon_pos(m, 2);
		
		if m == size(tx_exon_pos, 1)
			fprintf(fid, '%d],\n', left);
		else
			fprintf(fid, '%d, ', left);
		end
	end
	
	if size(tx_exon_pos, 1) == 0
		fprintf(fid, '1],\n');
	end
	
	cds = transcripts.CDS(tx_idx, :);
	if ~any(isnan(cds))
		cds(isnan(cds)) = 0;
		fprintf(fid, '\t"cds": [%d, %d],\n', cds(1), cds(2));
	end
	
	fprintf(fid, '\t"probe_pos": [ ');
	
	for m = 1:size(probe_pos, 1)
		if m == size(probe_pos, 1)
			fprintf(fid, '[%d, %d] ],\n', probe_pos(m, 1), probe_pos(m, 2));
		else
			fprintf(fid, '[%d, %d], ', probe_pos(m, 1), probe_pos(m, 2));
		end
    end
	
	if size(probe_pos, 1) == 0
		fprintf(fid, '],\n');
	end
	
	fprintf(fid, '\t"probe_accepted": [ ');
	for m = 1:size(probe_pos, 1)
		if m == size(probe_pos, 1)
			fprintf(fid, '%d ],\n', accepted(m));
		else
			fprintf(fid, '%d, ', accepted(m));
		end
	end
	
	if size(probe_pos, 1) == 0
		fprintf(fid, '],\n');
	end

	
	fprintf(fid, '\t"probe_seq": [ ');
	
	for m = 1:size(probe_pos, 1)
		if m == size(probe_pos, 1)
			fprintf(fid, '"%s"]\n', probe_seq{m});
		else
			fprintf(fid, '"%s", ', probe_seq{m});
		end
    end
	
	if size(probe_pos, 1) == 0
		fprintf(fid, ']\n');
	end
		
	if t == length(tx_indices)
		fprintf(fid, '}\n');
	else
		fprintf(fid, '},\n');
	end
end

fprintf(fid, '}\n');
fclose(fid);

