
% Author: Matti Annala <matti.annala@tut.fi>

function junctions = exon_connectivity(reads, read_len, varargin)
	
global organism;
genes = organism.Genes;
exons = organism.Exons;

max_mismatches = 2;
max_exon_skips = 2;
min_anchor = 10;

flank_len = read_len - min_anchor

drop_args = false(length(varargin), 1);
for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MaxMismatches')
		max_mismatches = varargin{k+1};
		drop_args(k:k+1) = true;
		continue;
	end
end
varargin = varargin(~drop_args);

S = length(reads.Raw);
prealloc = 10e6;

junctions = struct;
junctions.Exons = nan(prealloc, 2);
junctions.Sequence = cell(prealloc, 1);
junctions.Reads = zeros(prealloc, S);

exon_numbers = regexprep(exons.ID, '[^\d]+', '');
exon_numbers = str2double(exon_numbers);

E = 1;

for g = 1:length(genes.Name)
	gene_exons = find(exons.Gene == g);
	[~, order] = sort_nat(exons.ID(gene_exons));
	gene_exons = gene_exons(order);
	
	gene_exnum = exon_numbers(gene_exons);
	
	if length(gene_exnum) <= 1, continue, end
	
	run_ends = [ find(gene_exnum(1:end-1) ~= gene_exnum(2:end)); ...
		length(gene_exnum) ];
	run_lengths = diff([0; run_ends]);
	
	if length(run_lengths) <= 1, continue, end
		
	exon_groups = cell(1, length(run_lengths));
	
	pos = 1;
	for r = 1:length(run_lengths)
		exon_groups{r} = gene_exons(pos:pos+run_lengths(r)-1);
		pos = pos + run_lengths(r);
	end
	
	for rg = 2:length(exon_groups)
		for lg = max(1, rg - max_exon_skips - 1):rg-1
			%exon_pairs = nan(length(exon_groups{lg})*length(exon_groups{rg}),2);
			for x = exon_groups{lg}'
				for y = exon_groups{rg}'
					junctions.Exons(E, :) = [x y];
					
					left_seq = exons.Sequence{x};
					if length(left_seq) > flank_len
						left_seq = left_seq(end-flank_len+1:end);
					end
					right_seq = exons.Sequence{y};
					if length(right_seq) > flank_len
						right_seq = right_seq(1:flank_len);
					end
					
					junctions.Sequence{E} = [left_seq, right_seq];
					E = E + 1;
				end
			end
		end
	end
end

junctions = filter_struct(junctions, 1:E);


% Construct an index of the putative exon-exon junctions.
tmp_pool = FilePool;
fasta_tmp = tmp_pool.temp('index_fasta');
index_tmp = tmp_pool.temp('index');

write_seq_fasta(junctions, fasta_tmp);

color_option = '';
index_suffix = '';
%if color
%	color_option = '-C';
%	index_suffix = '_colorspace';
%end

index_name = [index_tmp index_suffix];

[status, out] = unix(sprintf('%s/tools/bowtie/bowtie-build %s %s %s', ...
	ppath, color_option, fasta_tmp, index_name));
if status ~= 0, error('Bowtie index construction failed: %s', out); end






for s = 1:S
	extracted = extract_reads(filter_query(reads, s));
	
	fprintf(1, 'Finding unaligned transcriptome reads for sample #%d...\n', s);
	al = align_reads(extracted, index_tmp, 'MaxMismatches', max_mismatches, ...
		'AllowAlignments', 1, 'Columns', 'target');
	
	for j = str2double(al.Target)'
		junctions.Reads(j, s) = junctions.Reads(j, s) + 1;
	end
end

