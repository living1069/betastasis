function dist = mirna_clusters()

global organism;
mirnas = organism.miRNA;
pre_mirnas = organism.pre_miRNA;

dist = Inf(length(pre_mirnas.Name));

for a = 1:length(pre_mirnas.Name)
	chr_a = pre_mirnas.Chromosome(a);
	strand_a = pre_mirnas.Strand(a);
	pos_a = pre_mirnas.Position(a);
	
	if isnan(chr_a) || any(isnan(pos_a)), continue, end
	
	for b = 1:length(pre_mirnas.Name)
		chr_b = pre_mirnas.Chromosome(b);
		strand_b = pre_mirnas.Strand(b);
		pos_b = pre_mirnas.Position(b);
		
		if isnan(chr_b) || any(isnan(pos_b)), continue, end
		if chr_a ~= chr_b, continue, end
		if strand_a ~= strand_b, continue, end   % FIXME: Check the strand?
		
		dist(a, b) = abs(pos_a - pos_b);
	end
end

for a = 1:length(pre_mirnas.Name), dist(a, a) = 0; end

vdist = squareform(dist);

figure;
vdist = vdist(vdist ~= Inf);
subplot(211); hist(vdist, 200);
vdist = vdist(vdist < 1e6);
subplot(212); hist(vdist, 200);
saveas(gcf, '~/premir_clusters.pdf');

conn = (dist < 1e5);
graph = biograph(conn);


[S, C] = conncomp(graph);

% Print the microRNA clusters.
for c = 1:S-1
	pre = find(C == c);
	if length(pre) <= 1, continue, end
		
	chr = pre_mirnas.Chromosome(pre(1));
	
	fprintf(1, 'Cluster #%d (chr%s):\n', c, organism.Chromosomes.Name{chr});
	for p = pre
		mir = pre_mirnas.Matures(p, 1:pre_mirnas.MatureCount(p));
		
		fprintf(1, '- %s: ', pre_mirnas.Name{p});
		fprintf(1, '%s', mirnas.Name{mir(1)});
		for m = mir(2:end), fprintf(1, ', %s', mirnas.Name{m}); end
		fprintf(1, '\n');
	end
end

