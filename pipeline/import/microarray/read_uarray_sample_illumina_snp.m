function [raw, probesets] = read_uarray_sample_beadstudio(sample_file)

global organism;

probe_id_type = '';

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
	
	tokens = regexp(line, 'Content\s+(.+)', 'tokens');
	if length(tokens) == 1
		token = tokens{1}; platform = token{1};
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
	if regexpi(platform, 'Cardio-Metabo_Chip')
		format = '%*s %s %*s %*s %*s %*s %f %s %d %f %f %*s %*s';
		data = textscan(fid, format, P, 'Delimiter', '\t');
		
		sample_ids = data{1};
		%gc_scores = data{2};
		chromosomes = chromosome_sym2num(data{3})';
		offsets = data{4};
		x_raw = data{5};
		y_raw = data{6};
		
	elseif regexpi(platform, 'Human660W-Quad')
		format = '%s %s %*s %*s %*s %*s %f %f %*s %f %*s';
		data = textscan(fid, format, P, 'Delimiter', '\t');
		
		probe_ids = data{1};
		sample_ids = data{2};
		x_raw = data{3};
		y_raw = data{4};
		
		if s == 1
			probes_fid = fopen( ...
			  '~/platforms/microarrays/illumina_660w_quad/ucsc_track_hg19.txt');
			probes = textscan(probes_fid, '%*s %s %d %*d %s %*s %*s %*s', ...
				'Delimiter', '\t', 'HeaderLines', 1);
			fclose(probes_fid);
			
			probe_id_to_idx = containers.Map(probes{3}, ...
				num2cell(1:length(probes{3})));
			
			valid_probes = probe_id_to_idx.isKey(probe_ids);
			idx = cell2mat(probe_id_to_idx.values(probe_ids(valid_probes)));
			chromosomes = chromosome_sym2num(probes{1}(idx));
			offsets = probes{2}(idx);
		end
		
		
		
	else
		error 'Unrecognized Illumina microarray platform.';
	end
	
	if any(~strcmp(sample_ids{1}, sample_ids))
		error 'Samples were mixed up.';
	end
	
	raw.Meta.Sample.Filename{s, 1} = sample_file;
	raw.Meta.Sample.ID{s, 1} = sample_ids{1};
	
	% Polar transformation to get the CNV intensity.
	raw.Mean(:, s) = sqrt(x_raw.^2 + y_raw.^2);
	
	if s == 1
		PS = sum(valid_probes);
		
		probesets = struct;
		probesets.Type = 'Copy number';
		probesets.Organism = 'Homo sapiens';
		probesets.GenomeVersion = 'GRCh37';
		probesets.Chromosome = chromosomes;
		probesets.Offset = offsets;
		probesets.ProbeCount = ones(PS, 1);
		probesets.Probes = find(valid_probes);
		
		% Order the probesets by chromosome and then by genomic offset.
		ordered_probecount = zeros(PS, 1);
		ordered_offset = zeros(PS, 1);
		ordered_probes = zeros(size(probesets.Probes));
		ordered_chromosome = zeros(PS, 1);

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
	end
	
	progress.update(s / S);
end

fclose(fid);


