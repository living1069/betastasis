
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
%    CREATE_GENE_PROBESETS(..., 'EveryTranscript', BOOL) lets the user specify
%    whether a probe is required to align to every transcript of a gene in order
%    to be included in its probeset. Enabled by default.
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
genome = organism.Genes;
transcriptome = organism.Transcripts;

include_genes = true(length(genome.Name), 1);
require_every_transcript = true;
min_probes_per_probeset = 4;
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

fprintf(1, 'Writing probe sequences to a temporary file...\n');
probes_fasta_tmp = ptemp;
write_probes_fasta(probes_fasta_tmp, probes);

fprintf(1, 'Aligning probes to transcripts using Bowtie...\n');
[alignments_tmp, ~] = bowtie_align(probes_fasta_tmp, 'transcripts', ...
	sprintf('-p4 -B1 -v%d -y --all --suppress 5,6,7,8', max_mismatches));

fprintf(1, 'Constructing gene expression probesets based on alignments...\n');

transcript_map = containers.Map(transcriptome.Name, ...
	num2cell(1:length(transcriptome.Name)));

probesets = struct;
probesets.Gene = genome.Name;
probesets.ProbeCount = zeros(length(genome.Name), 1);
probesets.Probes = zeros(length(genome.Name), 10);

% Read the Bowtie alignments into Matlab.
alignments_file = fopen(alignments_tmp);
data = textscan(alignments_file, '%d %*s %s %d');
fclose(alignments_file);

probe_indices = data{1};
transcripts = data{2};
offsets = data{3};

% Map the transcript names found in the Bowtie alignments file into transcript
% indices. It is faster to map them all in one go, rather than inside the loop.
transcript_indices = cell2mat(transcript_map.values(transcripts));

% Find sequential lines that indicate multiple alignments for one probe.
% These are identified by looking at the probe ID column.
run_ends = [ find(probe_indices(1:end-1) ~= probe_indices(2:end)); ...
	length(probe_indices) ];
run_lengths = diff([0; run_ends]);

if ~isempty(report_dir)
	report = struct;
	report.TxProbePosition = cell(length(organism.Transcripts.Name), 1);
	report.TxProbeAccepted = cell(length(organism.Transcripts.Name), 1);
	report.TxProbeSeq = cell(length(organism.Transcripts.Name), 1);
	
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
	probe_target_genes = unique(transcriptome.Gene(probe_target_transcripts));
	
	if require_every_transcript
		% For each potential target gene, check if the probe targets every
		% transcript variant of the gene.
		matched_genes = [];
		for g = 1:length(probe_target_genes)
			idx = probe_target_genes(g);
			if isempty(setdiff( ...
				genome.Transcripts(idx, 1:genome.TranscriptCount(idx)), ...
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
			probe_num = probesets.ProbeCount(idx);
			if probe_num < 40
				probesets.Probes(idx, probe_num + 1) = probe_indices(pos);
				probesets.ProbeCount(idx) = probe_num + 1;
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

if min_probes_per_probeset > 0
	fprintf(1, 'Filtering out probesets with less than %d probes...\n', ...
		min_probes_per_probeset);
	probesets.ProbeCount(probesets.ProbeCount < min_probes_per_probeset) = 0;
	probesets.Probes(probesets.ProbeCount < min_probes_per_probeset, :) = 0;
end

if any(~include_genes)
	probesets.Probes(~include_genes, :) = 0;
	probesets.ProbeCount(~include_genes) = 0;
end

probesets.Type = 'Gene expression';
probesets.Organism = organism.Name;
probesets.GenomeVersion = organism.GenomeVersion;

% Remove temporary files.
safe_delete(probes_fasta_tmp);
safe_delete(alignments_tmp);

% Finally we write a report about the generated probesets.
if ~isempty(report_dir)
	fprintf(1, 'Writing a probeset report to %s...\n', report_dir);
	
	[~, ~] = mkdir(report_dir);
	
	fid = fopen([report_dir '/genelist.json'], 'w');
	fprintf(fid, '{ "genes": [ ');
	first = true;
	
	for g = 1:length(organism.Genes.Name)
		if probesets.ProbeCount(g) == 0, continue, end
		plot_probe_reads(report, g, report_dir);
		if first 
			fprintf(fid, '"%s"', organism.Genes.Name{g});
			first = false;
		else
			fprintf(fid, ', "%s"', organism.Genes.Name{g});
		end
	end
	
	fprintf(fid, ' ] }\n');
	fclose(fid);
end








function [] = plot_probe_reads(report, gene_idx, report_dir)
	
global organism;

tx_indices = find(organism.Transcripts.Gene == gene_idx);

longest_tx = 0;
for t = 1:length(tx_indices)
	tx_seq = organism.Transcripts.Sequence{tx_indices(t)};
	longest_tx = max(longest_tx, length(tx_seq));
end

gene_name = organism.Genes.Name{gene_idx};
[~, ~] = mkdir([report_dir '/' lower(gene_name(1))]);

fid = fopen([report_dir '/' lower(gene_name(1)) '/' ...
	organism.Genes.Name{gene_idx} '.json'], 'w');
fprintf(fid, '{\n');

for t = 1:length(tx_indices)
	tx_idx = tx_indices(t);
	tx_seq = organism.Transcripts.Sequence{tx_idx};
	
	probe_pos = report.TxProbePosition{tx_idx};
	accepted = report.TxProbeAccepted{tx_idx};
	probe_seq = report.TxProbeSeq{tx_idx};
	
	fprintf(fid, '"%s": {\n', organism.Transcripts.Name{tx_idx});
	fprintf(fid, '\t"length": %d,\n', length(tx_seq));
	fprintf(fid, '\t"exons": [');
	
	tx_exons = find(organism.Exons.Transcript == tx_idx);
	for m = 1:length(tx_exons)
		left = organism.Exons.Position(tx_exons(m), 1);
		right = organism.Exons.Position(tx_exons(m), 2);
		
		if m == length(tx_exons)
			fprintf(fid, '%d],\n', left);
		else
			fprintf(fid, '%d, ', left);
		end
	end
	
	if isempty(tx_exons)
		fprintf(fid, '1],\n');
	end
	
	cds = organism.Transcripts.CDS(tx_idx, :);
	if ~any(isnan(cds))
		cds(isnan(cds)) = 0;
		fprintf(fid, '\t"cds": [%d, %d],\n', cds(1), cds(2));
	end
	
	read_count_dist = zeros(1, length(tx_seq));
	
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

