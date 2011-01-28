function [raw, probesets] = read_uarray_sample_beadstudio(sample_file)

global organism;

fid = fopen(sample_file);
while 1
	line = fgetl(fid);
	if line == -1, break, end
		
	tokens = regexp(line, 'Num SNPs\s+(\d+)', 'tokens');
	if length(tokens) == 1
		token = tokens{1};
		P = str2double(token{1});
	end
	
	tokens = regexp(line, 'Num Samples\s+(\d+)', 'tokens');
	if length(tokens) == 1
		token = tokens{1};
		S = str2double(token{1});
	end

	if strcmp(line, '[Data]'), break, end
end

header = fgetl(fid);

raw = struct;
raw.Mean = nan(P, S);
raw.Meta.Type = 'Microarray probe intensities';
raw.Meta.Platform = 'Illumina iSelect Cardio-MetaboChip';

progress = Progress;

for s = 1:S
	data = textscan(fid, '%*s %s %*s %*s %*s %*s %f %s %d %f %f %*s %*s', P, ...
		'Delimiter', '\t');
		
	sample_ids = data{1};
	gc_scores = data{2};
	chromosomes = data{3};
	offsets = data{4};
	x_raw = data{5};
	y_raw = data{6};
	
	if any(~strcmp(sample_ids{1}, sample_ids))
		error 'Samples were mixed up.';
	end
	
	raw.Meta.Sample.Filename{s, 1} = sample_file;
	raw.Meta.Sample.ID{s, 1} = sample_ids{1};
	
	% Polar transformation to get the CNV intensity.
	raw.Mean(:, s) = sqrt(x_raw.^2 + y_raw.^2);
	
	if s == 1
		probesets = struct;
		probesets.Type = 'Copy number';
		probesets.Organism = 'Homo sapiens';
		probesets.GenomeVersion = 'GRCh37';
		probesets.Chromosome = chromosome_sym2num(chromosomes)';
		probesets.Offset = offsets(1:P);
		probesets.ProbeCount = ones(P, 1);
		probesets.Probes = (1:P)';
		
		% Order the probesets by chromosome and then by genomic offset.
		ordered_probecount = zeros(P, 1);
		ordered_offset = zeros(P, 1);
		ordered_probes = zeros(size(probesets.Probes));
		ordered_chromosome = zeros(P, 1);

		pos = 1;
		for k = 1:length(organism.Chromosomes.Name)
			chr_psets = find(probesets.Chromosome == k);
			n = length(chr_psets);
			
			ordered_chromosome(pos:pos+n-1) = k;
			[ordered_offset(pos:pos+n-1), order] = ...
				sort(probesets.Offset(chr_psets));
			
			idx = chr_psets(order);
			ordered_probecount(pos:pos+n-1) = probesets.ProbeCount(idx);
			ordered_probes(pos:pos+n-1, :) = probesets.Probes(idx, :);
			
			pos = pos + n;
		end

		probesets.Chromosome = ordered_chromosome;
		probesets.Offset = ordered_offset;
		probesets.ProbeCount = ordered_probecount;
		probesets.Probes = ordered_probes;
		return;
	end
	
	progress.update(s / S);
end

fclose(fid);


