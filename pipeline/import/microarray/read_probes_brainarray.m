
% This function reads a set of probes and original probesets from a
% tab-delimited BrainArray probe description file.
% 
% Inputs:
%     filepath - Filesystem path to the probe description file.
%
% Outputs:
%     probes - Data structure that contains the probe positions and sequences.
%     probesets - Data structure that describes the BrainArray probesets.
%
% Author: Matti Annala <matti.annala@tut.fi>

function [probes, probesets] = read_probes_brainarray(filepath)

global organism;
transcriptome = organism.Transcripts;

[probes, orig_probesets] = read_probes_affy_tab(filepath);

fprintf(1, 'Translating original BrainArray probesets to common format...\n');

transcript_map = containers.Map(transcriptome.Name, ...
	num2cell(1:length(transcriptome.Name)));

probesets = struct( ...
	'Transcript', { transcriptome.Name }, ...
	'ProbeCount', zeros(length(transcriptome.Name), 1), ...
	'Probes', zeros(length(transcriptome.Name), 10), ...
	'GenomeVersion', organism.GenomeVersion);

for ps = 1:length(orig_probesets.ProbeCount)
	name = orig_probesets.Name{ps};
	if strcmp(name(end-2:end), '_at'), name = name(1:end-3); end
	if ~transcript_map.isKey(name)
		fprintf(1, 'Transcript name %s not found in RefSeq.\n', name);
		continue;
	end
	
	transcript = transcript_map(name);
	probe_num = orig_probesets.ProbeCount(ps);
	probesets.ProbeCount(transcript) = probe_num;
	probesets.Probes(transcript, 1:probe_num) = ...
		orig_probesets.Probes(ps, 1:probe_num);
end

