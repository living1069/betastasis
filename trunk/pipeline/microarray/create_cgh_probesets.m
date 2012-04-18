function probesets = create_cgh_probesets(probes)

global organism;

tmp = temporary('create_cgh_probesets');

probesets = struct;
probesets.Type = 'Copy number';
probesets.Organism = organism.Name;

fprintf(1, 'Aligning probes to genomic locations using Bowtie...\n');
write_seq_fasta(probes, [tmp 'probes.fa']);
alignments = bowtie_align(import_reads(tmp), 'genome', '-v0 -m1');

fprintf(1, 'Constructing copy number probesets based on alignments...\n');

al = all_alignments(alignments);
probe_indices = str2double(al.read);
chromosomes = al.target;
offsets = al.offset;

probeset_map = containers.Map;

progress = Progress;

for k = 1:length(probe_indices)
	probeset_name = [ chromosomes{k} '_' num2str(offsets(k)) ];
	if ~probeset_map.isKey(probeset_name)
		probeset_map(probeset_name) = {};
	end
	
	probeset_map(probeset_name) = cat(2, probeset_map(probeset_name), ...
		{ probe_indices(k) });
	progress.update(k / length(probe_indices));
end

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

