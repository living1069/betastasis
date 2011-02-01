
% CREATE_EXON_PROBESETS   Builds custom exon expression probesets from probes
% 
%    PROBESETS = CREATE_EXON_PROBESETS(PROBES, ...) constructs custom
%    exon expression probesets from the microarray probes described by the data
%    structure PROBES. The probesets are designed based on exon annotations
%    in the currently selected transcriptome build (see HELP ORGANISM).
%    
%    The probesets are built by aligning the probe sequences against the
%    reference exome. A probe is then assigned to the probesets of any exons
%    whose sequences it aligns against.
%
%    CREATE_GENE_PROBESETS(..., 'MaxMismatches', MAX_MM) allows MAX_MM 
%    nucleotide mismatches when aligning the probes against the exons.
%    Default is 0 nucleotide mismatches.
%
%    CREATE_GENE_PROBESETS(..., 'MinProbesetSize', NPROBES) discards probesets
%    with less than NPROBES probes. Default is NPROBES = 2.

% Author: Matti Annala <matti.annala@tut.fi>

function probesets = create_exon_probesets(probes, varargin)

global organism;
exons = organism.Exons;

min_probes_per_probeset = 2;
max_mismatches = 0;
allow_nonspecific_hits = false;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MaxMismatches')
		max_mismatches = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'MinProbesetSize')
		min_probes_per_probeset = varargin{k+1};
		continue;
	end
	
	%if strcmpi(varargin{k}, 'AllowNonspecific')
	%	allow_nonspecific_hits = varargin{k+1};
	%	continue;
	%end
	
	error('Unrecognized option "%s".', varargin{k});
end

fprintf(1, 'Writing probe sequences to a temporary file...\n');
probes_fasta_tmp = ptemp;
write_seq_fasta(probes, probes_fasta_tmp);

fprintf(1, 'Aligning probes to exons using Bowtie...\n');
[alignments_tmp, ~] = bowtie_align(probes_fasta_tmp, 'exons', ...
	sprintf('-v%d -m10 -k10 --suppress 5,6,7,8 ', max_mismatches));

fprintf(1, 'Constructing exon probesets based on alignments...\n');

probesets = struct;
probesets.ProbeCount = zeros(length(exons.ID), 1);
probesets.Probes = zeros(length(exons.ID), 4);

alignments_file = fopen(alignments_tmp);
data = textscan(alignments_file, '%d %*s %d %*d');
fclose(alignments_file);

probe_indices = data{1};
alignment_target_exons = data{2};
clear data;

for p = 1:length(probe_indices)
	ex = alignment_target_exons(p);
	probe_num = probesets.ProbeCount(ex);
	if probe_num < 20
		probesets.Probes(ex, probe_num + 1) = probe_indices(p);
		probesets.ProbeCount(ex) = probe_num + 1;
	end
end

if min_probes_per_probeset > 0
	fprintf(1, 'Filtering out probesets with less than %d probes...\n', ...
		min_probes_per_probeset);
	probesets.ProbeCount(probesets.ProbeCount < min_probes_per_probeset) = 0;
	probesets.Probes(probesets.ProbeCount < min_probes_per_probeset, :) = 0;
end

probesets.Type = 'Exon expression';
probesets.Organism = organism.Name;
%probesets.GenomeVersion = organism.GenomeVersion;

safe_delete(probes_fasta_tmp);
safe_delete(alignments_tmp);

