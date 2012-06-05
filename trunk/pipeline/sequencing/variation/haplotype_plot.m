
function haplo = haplotype_plot(variants)

global organism;
chromosomes = organism.Chromosomes;

S = length(variants.meta.sample_id);
window_size = 500e3;

% Filter out variants with too low genotype quality.
keep = all(variants.genotype_quality >= 50, 2);
variants = filter_rows(variants, keep);

% Filter out variants where there is no deviation from the reference genome.
keep = any(variants.genotype > 0, 2);
variants = filter_rows(variants, keep);

var = variants.rows;

haplo.chromosome = [];
haplo.position = [];
haplo.score = [];
haplo.total_variants = [];

window_first_var = [];
window = [1 window_size] - window_size;
chr = 1;

progress = Progress;

for k = 1:length(variants.rows.chromosome)
	if var.chromosome(k) ~= chr || var.position(k) > window(2)
		wvar = window_first_var:(k-1);
		if length(wvar) >= 10
			haplo.chromosome(end+1, 1) = chr;
			haplo.position(end+1, 1) = window(1) + round(window_size / 2);
			haplo.total_variants(end+1, 1) = length(wvar);
			haplo.score(end+1, 1) = mean(max( ...
				sum(variants.genotype(wvar, :) == 0, 2), ...
				sum(variants.genotype(wvar, :) > 0, 2)));
		end
		
		if var.chromosome(k) ~= chr
			progress.update(chr / length(chromosomes.Name));
			chr = var.chromosome(k);
			window = [1 window_size] - window_size;
		end
		
		while var.position(k) > window(2), window = window + window_size; end
		window_first_var = k;
	end
end

progress.update(1);

% Finally we normalize the conservation scores to the range [0,1].
haplo.score = haplo.score - ceil(S / 2);
haplo.score = haplo.score / (S - ceil(S / 2));

for c = 1:24
	figure;
	cidx = (haplo.chromosome == c);
	scatter(haplo.position(cidx), haplo.score(cidx));
	ylim([0 1]);
	chromosomeplot('/data/csb/organisms/homo_sapiens/cytoBand.txt', c, ...
		'addtoplot', gca, 'Orientation', 2);
	saveas(gcf, sprintf('~/chr%s_haplo.pdf', chromosomes.Name{c}));
end




% Draw a plot that shows all conserved regions across the whole genome.
cons = struct;
cons.Chromosome = [];
cons.CNVType = [];
cons.Start = [];
cons.End = [];

for c = 1:24
	conserved = (haplo.score > 0.95 & haplo.chromosome == c);
	conserved = medfilt1(double(conserved), 5);
	
	run_ends = [ find(conserved(1:end-1) ~= conserved(2:end)); ...
		length(conserved) ];
	run_lengths = diff([0; run_ends]);
	
	pos = 1;
	for r = 1:length(run_lengths)
		if conserved(pos)
			cons.Chromosome(end+1, 1) = c;
			cons.CNVType(end+1, 1) = 2;
			cons.Start(end+1, 1) = haplo.position(pos);
			cons.End(end+1, 1) = haplo.position(pos+run_lengths(r)-1);
		end
		pos = pos + run_lengths(r);
	end
end

figure; chromosomeplot('/data/csb/organisms/homo_sapiens/cytoBand.txt', ...
	'CNV', cons);
saveas(gcf, '~/all_chr_haplo.pdf');





% Print the genes found within the conserved regions.
genes = organism.Genes;
cons_genes = {};
for c = find(cons.Chromosome == 23)'
	cons_genes = [cons_genes; genes.Name( ...
		genes.Chromosome == cons.Chromosome(c) & ...
		genes.Position(:, 1) > cons.Start(c) & ...
		genes.Position(:, 2) < cons.End(c))];
end

cons_genes = unique(cons_genes);
fprintf('CONSERVED GENES:\n');
fprintf('%s\n', cons_genes{:});


