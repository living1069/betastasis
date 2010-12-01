function [] = select_organism(name, version)

name(1) = upper(name(1));

orgpath = [ppath '/organisms/' flatten_str(name) '/' flatten_str(version)];
if exist(orgpath) ~= 7
	fprintf(1, 'Could not find the specified organism.\n');
	return;
end

global organism;
organism = struct;
organism.Name = name;
organism.Version = version;
organism.Chromosomes = Chromosomes;
organism.Genes = Genes;
organism.Transcripts = Transcripts;
organism.miRNA = MiRNA;
organism.pre_miRNA = Pre_miRNA;

