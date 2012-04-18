
% Author: Matti Annala <matti.annala@tut.fi>

function [] = write_bigwig_track(val, pos, track_file, varargin)

global organism;
chromosomes = organism.Chromosomes;

%tmp = temporary('write_bigwig_track');

step = pos.position(2) - pos.position(1);
chr_starts = [0 find(pos.position(1:end-1) > pos.position(2:end))'] + 1;
chr_ends = [chr_starts(2:end)-1 length(pos.position)];

chr_ends = chr_ends - 1;    % FIXME: To make wigToBigWig work.

fid = fopen(track_file, 'W');
for c = 1:length(chr_starts)
	fprintf(fid, 'fixedStep chrom=chr%s start=%d step=%d\n', ...
		chromosomes.Name{c}, pos.position(chr_starts(c)), step);
	fprintf(fid, '%f\n', val(chr_starts(c):chr_ends(c)));
end

fclose(fid);

