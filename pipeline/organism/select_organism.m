function [] = select_organism(organism_name)

orgpath = [ppath '/organisms/' flatten_str(organism_name)];

fprintf(1, 'Reading organism data into memory...\n');

global organism;
organism = load([orgpath '/default']);
organism = organism.organism;

