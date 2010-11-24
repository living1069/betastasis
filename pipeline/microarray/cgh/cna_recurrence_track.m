function [] = cna_recurrence_track(samples, refs, cgh_probesets, track_file, ...
	varargin)

global organism;

cna = cna_from_cgh(samples, refs, cgh_probesets, varargin{:});

S = size(cna, 2);
amp_recurrence = sum(cna > 0, 2) / S;
del_recurrence = sum(cna < 0, 2) / S;

cna_sign = ((amp_recurrence >= del_recurrence) - 0.5) * 2;
recurrence = cna_sign .* max(amp_recurrence, del_recurrence);

fprintf(1, 'Writing CNA recurrence profile as an IGV track...\n');

fid = fopen(track_file, 'W');
fprintf(fid, 'Chromosome\tStart\tEnd\tFeature\tCNA_Recurrence\n');

% Find all probes that target the chromosome we're interested in.
for chr = 1:24
	idx = find(cgh_probesets.Chromosome == chr);
	N = length(idx);
	
	offsets = cgh_probesets.Offset(idx);
	chr_recurrence = recurrence(idx);
	
	borders = round(mean([offsets(1:N-1) offsets(2:N)], 2));
	for p = 2:N-1
		fprintf(fid, '%d\t%d\t%d\t-\t%f\n', chr, borders(p-1)+1, borders(p), ...
			chr_recurrence(p));
	end
end

fclose(fid);

