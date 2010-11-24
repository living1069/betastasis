function probesets = create_cgh_probesets(probes)

global organism;

probesets = struct;
probesets.Type = 'Copy number';
probesets.Organism = organism.Name;
probesets.GenomeVersion = organism.GenomeVersion;

fprintf(1, 'Writing probe sequences to a temporary file...\n');
probes_fasta_tmp = ptemp();
write_seq_fasta(probes, probes_fasta_tmp);

fprintf(1, 'Aligning probes to genomic locations using Bowtie...\n');
[alignments_tmp, ~] = bowtie_align(probes_fasta_tmp, 'genome', ...
	'-v0 -m1 --suppress 5,6,7,8');

fprintf(1, 'Constructing copy number probesets based on alignments...\n');

alignments_file = fopen(alignments_tmp);
data = textscan(alignments_file, '%d %*s %s %d');
fclose(alignments_file);

probe_indices = data{1};
chromosomes = data{2};
offsets = data{3};
clear data;

probeset_map = containers.Map;

progress = 0;
fprintf(1, 'Progress: 00%%');

for k = 1:length(probe_indices)
	probeset_name = [ chromosomes{k} '_' num2str(offsets(k)) ];
	if ~probeset_map.isKey(probeset_name)
		probeset_map(probeset_name) = {};
	end
	
	probeset_map(probeset_name) = cat(2, probeset_map(probeset_name), ...
		{ probe_indices(k) });
	
	if floor(k / length(probe_indices) * 100) > progress
		progress = floor(k / length(probe_indices) * 100);
		fprintf(1, '\b\b\b%02d%%', progress);
	end
end

fprintf(1, '\n');

probeset_names = probeset_map.keys();

probesets.Chromosome = zeros(length(probeset_names), 1);
probesets.Offset = zeros(length(probeset_names), 1);
probesets.ProbeCount = zeros(length(probeset_names), 1);
probesets.Probes = zeros(length(probeset_names), 1);

for k = 1:length(probeset_names)
	name = probeset_names{k};
	tmp = find(name == '_');
	chr = chromosome_sym2num(name(1:(tmp(1)-1)));
	offset = str2num(name((tmp(1)+1):end));
	
	indices = probeset_map(name);
	
	probesets.Chromosome(k) = chr;
	probesets.Offset(k) = offset;
	probesets.ProbeCount(k) = length(indices);
	probesets.Probes(k, 1:probesets.ProbeCount(k)) = cell2mat(indices);
end

% Order the probesets by chromosome and then by genomic offset.
ordered_probecount = zeros(length(probeset_names), 1);
ordered_offset = zeros(length(probeset_names), 1);
ordered_probes = zeros(size(probesets.Probes));
ordered_chromosome = zeros(length(probeset_names), 1);

pos = 1;
for k = 1:length(organism.Chromosomes.Name)
	chr_psets = find(probesets.Chromosome == k);
	n = length(chr_psets);
	
	ordered_chromosome(pos:pos+n-1) = k;
	[ordered_offset(pos:pos+n-1), order] = sort(probesets.Offset(chr_psets));
	
	idx = chr_psets(order);
	ordered_probecount(pos:pos+n-1) = probesets.ProbeCount(idx);
	ordered_probes(pos:pos+n-1, :) = probesets.Probes(idx, :);
	
	pos = pos + n;
end

probesets.Chromosome = ordered_chromosome;
probesets.Offset = ordered_offset;
probesets.ProbeCount = ordered_probecount;
probesets.Probes = ordered_probes;

safe_delete(probes_fasta_tmp);
safe_delete(alignments_tmp);

