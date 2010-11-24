function probesets = create_mirna_probesets(probes)

global organism;
mirnas = organism.miRNA;

fprintf(1, 'Writing probe sequences to a temporary file...\n');
probes_fasta_tmp = ptemp;
write_seq_fasta(probes, probes_fasta_tmp);

alignments_tmp = ptemp;

fprintf(1, 'Aligning probes to transcripts using Bowtie...\n');
[alignments_tmp, ~] = bowtie_align(probes_fasta_tmp, 'mirnas', ...
	'-v0 -m1 --suppress 2,4,5,6,7,8');

fprintf(1, 'Constructing probesets based on Bowtie alignments...\n');

probesets = struct( ...
	'miRNA', { mirnas.Name }, ...
	'ProbeCount', zeros(length(mirnas.Name), 1), ...
	'Probes', zeros(length(mirnas.Name), 2));

alignments_file = fopen(alignments_tmp);
data = textscan(alignments_file, '%d %s');
fclose(alignments_file);

probe_indices = data{1};
mirna_indices = mirna_idx(data{2});
clear data;

for k = 1:length(probe_indices)
	idx = mirna_indices(k);
	probe_num = probesets.ProbeCount(idx);
	probesets.Probes(idx, probe_num + 1) = probe_indices(k);
	probesets.ProbeCount(idx) = probe_num + 1;
end

% Agilent miRNA arrays usually have two different types of probes for each
% type of microRNA.
%for k = 1:length(probesets.ProbeCount)
%	if probesets.ProbeCount(k) == 0, continue, end
%
%	seqs = probes.Sequence(probesets.Probes(k, 1:probesets.ProbeCount(k)));
%	useqs = unique(seqs);
%	
%	useq_lens = zeros(length(useqs), 1);
%	for m = 1:length(useqs), useq_lens(m) = length(useqs{m}); end
%		
%	[~, longest_useq] = max(useq_lens);
%	longest_seq = useqs(longest_useq(1));
%	
%	best_probes = find(strcmp(longest_seq, seqs));
%	probesets.ProbeCount(k) = length(best_probes);
%	probesets.Probes(k, 1:probesets.ProbeCount(k)) = ...
%		probesets.Probes(k, best_probes);
%end

%fprintf(1, 'Filtering out probesets with less than three probes...\n');
%probesets.ProbeCount(probesets.ProbeCount < 3) = 0;
%probesets.Probes(probesets.ProbeCount < 3, :) = 0;

probesets.Organism = organism.Name;
probesets.miRNAVersion = organism.miRNAVersion;
probesets.Type = 'miRNA expression';

system(['rm ' probes_fasta_tmp]);
system(['rm ' alignments_tmp]);

