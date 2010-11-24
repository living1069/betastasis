
% This function constructs exon expression probesets for the current organism,
% using a set of complementary DNA probes specified by the user.
%
% The algorithm is not very smart, and simply places a probe in an exon's
% probeset if the probe's target sequence or its reverse complement is found
% in the exon. The biggest issue with this is that some probes will be placed
% in multiple probesets.
%
% Inputs:
%     probes - Data structure describing a set of complementary DNA probes.
%     probe_blacklist - A list of probes that should be ignored and not placed
%	      in any probesets. The probes should be listed by their index.
% 
% Outputs:
%     probesets - Data structure describing a set of probesets, each targeting
%         a specific exon of the currently selected organism.
%
% Author: Matti Annala <matti.annala@tut.fi>

function probesets = create_exon_probesets(probes)

global organism;
exons = organism.Exons;

fprintf(1, 'Writing probe sequences to a temporary file...\n');
probes_fasta_tmp = ptemp();
write_seq_fasta(probes, probes_fasta_tmp);

fprintf(1, 'Aligning probes to transcripts using Bowtie...\n');
[alignments_tmp, ~] = bowtie_align(probes_fasta_tmp, 'exons', ...
	'-v0 -m10 --suppress 5,6,7,8 ');

fprintf(1, 'Constructing exon probesets based on alignments...\n');

probesets = struct( ...
	'ProbeCount', zeros(length(exons.Transcript), 1), ...
	'Probes', zeros(length(exons.Transcript), 4));

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

fprintf(1, 'Filtering out probesets with less than three probes...\n');
probesets.ProbeCount(probesets.ProbeCount < 3) = 0;
probesets.Probes(probesets.ProbeCount < 3, :) = 0;

probesets.Type = 'Exon expression';
probesets.Organism = organism.Name;
probesets.GenomeVersion = organism.GenomeVersion;

system(['rm ' probes_fasta_tmp]);
system(['rm ' alignments_tmp]);

