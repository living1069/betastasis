
% This function constructs transcript expression probesets for the current
% organism, using a set of complementary DNA probes specified by the user.
%
% The algorithm is not very smart, and simply places a probe in a transcript's
% probeset if the probe's target sequence or its reverse complement is found
% in the transcript. The biggest issue with this is that some probes will be
% placed in multiple probesets.
%
% Inputs:
%     probes - Data structure describing a set of complementary DNA probes.
%     probe_blacklist - A list of probes that should be ignored and not placed
%	      in any probesets. The probes should be listed by their index.
% 
% Outputs:
%     probesets - Data structure describing a set of probesets, each targeting
%         a specific transcript of the currently selected organism.
%
% Author: Matti Annala <matti.annala@tut.fi>

function probesets = create_transcript_probesets(probes)

global organism;
transcriptome = organism.Transcripts;

fprintf(1, 'Writing probe sequences to a temporary file...\n');
probes_fasta_tmp = ptemp;
write_seq_fasta(probes, probes_fasta_tmp);

fprintf(1, 'Aligning probes to transcripts using Bowtie...\n');
[alignments_tmp, ~] = bowtie_align(probes_fasta_tmp, 'transcripts', ...
	'-v0 -m10 --suppress 5,6,7,8');

fprintf(1, 'Constructing probesets based on Bowtie alignments...\n');

transcript_map = containers.Map(transcriptome.Name, ...
	num2cell(1:length(transcriptome.Name)));

probesets = struct( ...
	'Transcript', { transcriptome.Name }, ...
	'ProbeCount', zeros(length(transcriptome.Name), 1), ...
	'Probes', zeros(length(transcriptome.Name), 10));

alignments_file = fopen(alignments_tmp);
data = textscan(alignments_file, '%d %s %s %*d');
fclose(alignments_file);

probe_indices = data{1};
transcripts = data{3};
clear data;

% Map the transcript names found in the Bowtie alignments file into transcript
% indices. It is faster to map them all in one go, rather than inside the loop.
transcript_indices = cell2mat(transcript_map.values(transcripts));

% Find sequential lines that indicate multiple alignments for one probe.
% These are identified by looking at the probe ID column.
run_ends = [ find(probe_indices(1:end-1) ~= probe_indices(2:end)); ...
	length(probe_indices) ];
run_lengths = diff([0; run_ends]);

pos = 1;
for r = 1:length(run_lengths)
	% Determine transcript IDs for the transcripts that the probe matches with.
	probe_target_transcripts = transcript_indices(pos:pos+run_lengths(r)-1);
	
	% Determine a set of genes so that at least one of the target transcripts
	% belongs to each gene in the set.
	probe_target_genes = unique(transcriptome.Gene(probe_target_transcripts));
	
	% If the probe matches transcripts from multiple different genes, don't
	% insert the probe into any probeset, as it is uninformative.
	if length(probe_target_genes) == 1
		for k = 1:length(probe_target_transcripts)
			idx = probe_target_transcripts(k);
			probe_num = probesets.ProbeCount(idx);
			if probe_num < 20
				probesets.Probes(idx, probe_num + 1) = probe_indices(pos);
				probesets.ProbeCount(idx) = probe_num + 1;
			end
		end
	end
	
	pos = pos + run_lengths(r);
end

fprintf(1, 'Filtering out probesets with less than three probes...\n');
probesets.ProbeCount(probesets.ProbeCount < 3) = 0;
probesets.Probes(probesets.ProbeCount < 3, :) = 0;


for p = 1:length(probe_indices)
	ps = transcript_map(transcripts{p});
	probe_num = probesets.ProbeCount(ps);
	probesets.Probes(ps, probe_num + 1) = probe_indices(p);
	probesets.ProbeCount(ps) = probe_num + 1;
end

probesets.Organism = organism.Name;
probesets.Version = organism.Version;
probesets.Type = 'Transcript expression';

system(['rm ' probes_fasta_tmp]);
system(['rm ' alignments_tmp]);

