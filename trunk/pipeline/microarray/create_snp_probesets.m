function probesets = create_snp_probesets(probes)

global organism;

probesets = struct;
probesets.Type = 'SNP';
probesets.Organism = organism.Name;
probesets.GenomeVersion = 'hg19';

fprintf(1, 'Writing probe sequences to a temporary file...\n');
probes_fasta_tmp = [ptemp '.fa'];
write_seq_fasta(probes, probes_fasta_tmp);

fprintf(1, 'Aligning probes to genomic locations using Bowtie...\n');
al = align_reads(probes_fasta_tmp, 'genome', 'MaxMismatches', 1, ...
	'AllowAlignments', 1, 'Columns', 'read,target,offset,mismatches');

safe_delete(probes_fasta_tmp);

fprintf(1, 'Constructing SNP probesets based on alignments...\n');

progress = Progress;

probe_indices = al.ReadID;
chromosomes = al.Target;
offsets = al.Offset;
mismatches = al.Mismatches;

probeset_names = cell(length(probe_indices), 1);
for k = 1:length(probe_indices)
	probeset_names{k} = sprintf('%s_%d', chromosomes{k}, offsets(k));
	progress.update(k / length(probe_indices) * 0.2);
end

[probeset_names, ~, groups] = unique(probeset_names);

N = length(probeset_names) * 4;

probesets.Chromosome = zeros(N, 1);
probesets.Offset = zeros(N, 1);
probesets.RefAllele = repmat(' ', N, 1);
probesets.MutAllele = repmat(' ', N, 1);
probesets.RefProbes = zeros(N, 1);
probesets.MutProbes = zeros(N, 1);

valid_count = 0;

for k = 1:length(probeset_names)
	name = probeset_names{k};
	probes = find(groups == k);
	
	ref_probes = probes(strcmp('', mismatches(probes)));
	mut_probes = setdiff(probes, ref_probes);
	
	if isempty(ref_probes) || isempty(mut_probes), continue, end
	
	tmp = find(name == '_');
	chr = chromosome_sym2num(name(1:(tmp(1)-1)));
	offset = str2double(name((tmp(1)+1):end));
	
	[mut_signatures, ~, classes] = unique(mismatches(mut_probes));
	for s = 1:length(mut_signatures)
		sig_probes = mut_probes(classes == s);
		
		tokens = regexp(mut_signatures{s}, '(\d+):(.>.)', 'tokens');
		if length(tokens) ~= 1, continue, end
			
		token = tokens{1};
		sig_offset = str2double(token{1});
		sig = token{2};
		
		valid_count = valid_count + 1;
		probesets.Chromosome(valid_count) = chr;
		probesets.Offset(valid_count) = offset + sig_offset;
		probesets.RefAllele(valid_count) = sig(1);
		probesets.MutAllele(valid_count) = sig(3);
		probesets.RefProbes(valid_count, 1:length(ref_probes)) = ref_probes';
		probesets.MutProbes(valid_count, 1:length(mut_probes)) = mut_probes';
	end
	
	progress.update(0.2 + k / length(probeset_names) * 0.8);
end

probesets.Chromosome = probesets.Chromosome(1:valid_count);
probesets.Offset = probesets.Offset(1:valid_count);
probesets.RefAllele = probesets.RefAllele(1:valid_count);
probesets.MutAllele = probesets.MutAllele(1:valid_count);
probesets.RefProbes = probesets.RefProbes(1:valid_count, :);
probesets.MutProbes = probesets.MutProbes(1:valid_count, :);

fprintf(1, 'Ordering SNP probesets by chromosome and genomic offset.\n');
order = [];
for k = 1:length(organism.Chromosomes.Name)
	chr_psets = find(probesets.Chromosome == k);
	[~, chr_order] = sort(probesets.Offset(chr_psets));
	
	order = [order; chr_psets(chr_order)];
end

probesets.Chromosome = probesets.Chromosome(order);
probesets.Offset = probesets.Offset(order);
probesets.RefAllele = probesets.RefAllele(order);
probesets.MutAllele = probesets.MutAllele(order);
probesets.RefProbes = probesets.RefProbes(order, :);
probesets.MutProbes = probesets.MutProbes(order, :);

