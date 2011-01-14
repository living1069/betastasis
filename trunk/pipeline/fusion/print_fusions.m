function [] = print_fusions(fusions, varargin)

global organism;

min_reads = 0;
blacklist = {};

% Organism specific genes that are blacklisted by default.
if strcmpi(organism.Name, 'Homo sapiens')
	blacklist = {'RPPH1', 'LOC100008588', 'RN18S1', 'SNORD', 'SNORA'};
end

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MinReads')
		min_reads = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Blacklist')
		blacklist = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

[~, sort_indices] = sort(fusions.ReadCount, 1, 'descend');

fprintf(1, 'Potential fusion transcripts (ordered by prevalence):\n');
for k = 1:length(sort_indices)
	idx = sort_indices(k);
	
	sequences = {};
	if isfield(fusions, 'ReadSequences')
		sequences = fusions.ReadSequences(idx, 1:fusions.ReadCount(idx));
	end
	
	if fusions.ReadCount(idx) < min_reads, continue, end
	
	print_exon_pair(fusions.Exons(idx, 1), fusions.Exons(idx, 2), ...
		fusions.ReadCount(idx), sequences, blacklist);
end

return;



function [] = print_exon_pair(left_exon, right_exon, read_count, sequences, ...
	blacklist)

global organism;
genes = organism.Genes;
exons = organism.Exons;

left_gene = exons.Gene(left_exon);
right_gene = exons.Gene(right_exon);

for k = 1:length(blacklist)
	if regexp(genes.Name{left_gene}, blacklist{k}), return, end
	if regexp(genes.Name{right_gene}, blacklist{k}), return, end
end

fprintf(1, '- fusion of %s[%s] and %s[%s]:\n', genes.Name{left_gene}, ...
	exons.ID{left_exon}, genes.Name{right_gene}, exons.ID{right_exon});
fprintf(1, '  * supported by %d reads:\n', read_count);

for k = 1:length(sequences)
	fprintf(1, '    * %s\n', sequences{k});
end

return;
