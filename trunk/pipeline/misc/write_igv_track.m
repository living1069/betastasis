
% Author: Matti Annala <matti.annala@tut.fi>

function [] = write_igv_track(val, samples, pos, track_file, varargin)

global organism;
chromosomes = organism.Chromosomes;

range = [];

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Range')
		range = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

custom_format = false;

fid = fopen(track_file, 'W');
fprintf(fid, ...
	'#track maxHeightPixels=500:400:300 graphType=points\n');
if custom_format
	fprintf(fid, '#columns chr=1 start=2 data=3\n');
else
	fprintf(fid, 'Chromosome\tStart\tEnd\tFeature');
end

S = size(val, 2);
for s = 1:S, fprintf(fid, '\t%s', samples{s}); end
fprintf(fid, '\n');

% Write the logratio tracks into a file.
if custom_format
	for k = 1:size(val, 1)
		fprintf(fid, '%s\t%d', chromosomes.Name{pos.chromosome(k)}, ...
			pos.position(k));
		fprintf(fid, '\t%f', val(k, :));
		fprintf(fid, '\n');
	end
else
	for k = 1:size(val, 1)
		fprintf(fid, '%s\t%d\t%d\t-', chromosomes.Name{pos.chromosome(k)}, ...
			pos.position(k), pos.position(k));
		fprintf(fid, '\t%f', val(k, :));
		fprintf(fid, '\n');
	end
end


fclose(fid);

