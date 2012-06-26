function [] = write_bed(features, bed_file)

global organism;
chromosomes = organism.Chromosomes;

fid = fopen(bed_file, 'W');
for k = 1:length(features.chromosome)
	fprintf(fid, '%s\t%d\t%d\t%s\n', ...
		chromosomes.Name{features.chromosome(k)}, ...
		features.position(k, 1), features.position(k, 2), ...
		features.name{k});
end
fclose(fid);

