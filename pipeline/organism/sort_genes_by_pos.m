function sorted = sort_genes_by_pos(unsorted)

global organism;
genes = organism.Genes;

tmp = [genes.Chromosome(unsorted), genes.Position(unsorted), unsorted(:)];
tmp = sortrows(tmp, [1 2]);
sorted = tmp(:, 3);

