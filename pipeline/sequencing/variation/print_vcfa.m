
function [] = print_vcfa(variants, varargin)

global organism;
chromosomes = organism.Chromosomes;
genes = organism.Genes;

min_reads = 10;
discard_silent = false;
negative_samples = [];
ranking = 'position';
output_file = 'variants.vcfa';

for k = 1:2:length(varargin)
	if regexpi(varargin{k}, 'min.*reads')
		min_reads = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'discard.*silent')
		discard_silent = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'neg.*sample')
		negative_samples = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'output.*(file)?')
		output_file = varargin{k+1};
		continue;
	end
	
	if regexpi(varargin{k}, 'order|ranking')
		ranking = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end



V = length(variants.Alleles);
S = length(variants.Samples);

keep = true(V, 1);

% Filter out variants that show up in negative control samples.
if ~isempty(negative_samples)
	keep = keep & ~any(variants.Genotype(:, negative_samples) > 0, 2);
end

% Filter out genotypes that are based on less than the minimum amount of reads.
if ~isempty(min_reads)
	variants.Genotype(variants.TotalReads < min_reads) = NaN;
	keep = keep & ~all(isnan(variants.Genotype), 2);
end

% Filter out silent variants.
if discard_silent
	for v = find(keep)'
		missense = strfind(variants.Consequences{v}, 'Missense');
		if all(cellfun(@isempty, missense)), keep(v) = false; end
	end
end

V = sum(keep);



if regexpi(ranking, 'position')
	% No need to do anything, variants should already be ranked by position.
	order = 1:V;
elseif regexpi(ranking, 'frequency')
	% Sort the variants based on their frequency in the valid samples.
	score = nansum(variants.Genotype, 2);
	[~, order] = sort(score, 'descend');
	
else
	error('Unrecognized ranking method "%s".', ranking);
end

order = order(keep(order));





variants.Chromosome = variants.Chromosome(order);
variants.Position = variants.Position(order);
variants.Alleles = variants.Alleles(order);
variants.OverlapGenes = variants.OverlapGenes(order);
variants.ReferenceSeq = variants.ReferenceSeq(order);
variants.Consequences = variants.Consequences(order);
variants.Genotype = variants.Genotype(order, :);
variants.TotalReads = variants.TotalReads(order, :);





fid = fopen(output_file, 'W');
fprintf(fid, ['OVERLAP_GENES\tCONSEQUENCES\tREF_SEQUENCE\tCHROM\tPOS\t' ...
	'REF>ALT\tNUM_MUTATED\tNUM_HOMOZYGOUSLY_MUTATED']);
for s = 1:S, fprintf(fid, '\t%s', variants.Samples{s}); end
fprintf(fid, '\n');


for y = 1:length(variants.Alleles)
	% Write column with names of overlapped genes.
	fprintf_cell(fid, genes.Name(variants.OverlapGenes{y}), ', ');
	fprintf(fid, '\t');
	
	% Write column with SNV consequences.
	fprintf_cell(fid, variants.Consequences{y}, '. ');
	
	fprintf(fid, '.\t%s\tchr%s\t%d\t%s', variants.ReferenceSeq{y}, ...
		chromosomes.Name{variants.Chromosome(y)}, variants.Position(y), ...
		variants.Alleles{y});
		
	fprintf(fid, '\t%d\t%d', sum(variants.Genotype(y, :) > 0, 2), ...
		sum(variants.Genotype(y, :) > 1, 2));
	
	genotype_symbols = {'0/0', '1/0', '1/1'};
	for s = 1:S
		if ~isnan(variants.Genotype(y, s))
			fprintf(fid, '\t%s (%d)', ...
				genotype_symbols{variants.Genotype(y, s) + 1}, ...
				variants.TotalReads(y, s));
		else
			fprintf(fid, '\t');
		end
	end
	fprintf(fid, '\n');
end
fclose(fid);

