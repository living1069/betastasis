
% Author: Matti Annala <matti.annala@tut.fi>

function probesets = read_probesets_brainarray(filepath, affy_probes)

global organism;
transcripts = organism.Transcripts;

[probes, orig_probesets] = read_probes_affy_tab(filepath);

transcript_map = containers.Map(transcripts.Name, ...
	num2cell(1:length(transcripts.Name)));

probesets = struct;
probesets.Transcript = transcripts.Name;
probesets.ProbeCount = zeros(length(transcripts.Name), 1);
probesets.Probes = zeros(length(transcripts.Name), 10);
probesets.Organism = [organism.Name ' / ' organism.Version];

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

% Brainarray files are missing some discarded probes. To tie the probes with the
% original full list of probes, we need to map them based on positions.
max_rows = max(max(affy_probes.YPos), max(probes.YPos));
max_cols = max(max(affy_probes.XPos), max(probes.XPos));

affy_probe_idx = zeros(max_rows, max_cols);
affy_probe_idx(sub2ind(size(affy_probe_idx), affy_probes.YPos, ...
	affy_probes.XPos)) = 1:length(affy_probes.Sequence);

ba_to_affy = affy_probe_idx( ...
	sub2ind(size(affy_probe_idx), probes.YPos, probes.XPos));
	
% We want zeroes to map back to zeroes, so we augment the translation table.
ba_to_affy = [0; ba_to_affy];
probesets.Probes(:) = ba_to_affy(probesets.Probes(:) + 1);

