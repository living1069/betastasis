
% Author: Matti Annala <matti.annala@tut.fi>

function [] = write_bigwig_track(coverage, track_file)

global organism;
chromosomes = organism.Chromosomes;

rows = coverage.rows;
step = rows.position(2) - rows.position(1);
for c = 1:24
	chr_starts(c) = find(rows.chromosome == c, 1);
	chr_ends(c) = find(rows.chromosome == c, 1, 'last');
end

fid = fopen(track_file, 'W');
track_file
fid
for c = 1:24
	fprintf(fid, 'fixedStep chrom=chr%s start=%d step=%d\n', ...
		chromosomes.Name{c}, round(rows.position(chr_starts(c))), step);
	fprintf(fid, '%f\n', coverage.mean(chr_starts(c):chr_ends(c), 1));
end
fclose(fid);

