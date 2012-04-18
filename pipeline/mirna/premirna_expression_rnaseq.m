
% Author: Matti Annala <matti.annala@tut.fi>

function pre_expr = premirna_expression_rnaseq(reads, varargin)

global organism;
pre_mirnas = organism.pre_miRNA;

trim_len = nan;
max_mismatches = 2;

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'max.*mismatch')
		max_mismatches = varargin{k+1};
		continue;
	end
	
	if rx(varargin{k}, 'trim.*len')
		trim_len = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

S = length(reads.url);

pre_mirna_name_to_idx = containers.Map(pre_mirnas.Name, ...
	num2cell(1:length(pre_mirnas.Name)));

pre_expr.mean = zeros(length(pre_mirnas.Name), S);

pre_expr.meta = reads.meta;
pre_expr.meta.type = 'pre-miRNA expression';
pre_expr.meta.organism = organism.Name;
pre_expr.meta.mismatches_allowed = ones(1, S) * max_mismatches;
pre_expr.meta.trim_length = ones(1, S) * trim_len;
pre_expr.meta.total_reads = nan(1, S);
pre_expr.meta.total_aligned_reads = nan(1, S);

for s = 1:S
	fprintf('Calculating pre-microRNA expression levels for sample %s:\n', ...
		reads.meta.sample_id{s});
	
	if trim_len > 0
		trimmed = trim_reads(filter(reads, s), trim_len);
	else
		trimmed = filter(reads, s);
	end
	
	fprintf('-> Aligning reads to miRBase pre-miRNAs using Bowtie...\n');
	alignments = bowtie_align(trimmed, 'pre_mirnas', ...
		sprintf('-v%d -m1', max_mismatches));
	
	pre_expr.meta.total_reads(s) = alignments.total_reads;
	
	fprintf(1, '-> Summarizing pre-miRNA expression levels...\n');
	
	for al = iterate_alignments(alignments)
		pre_idx = cell2mat(pre_mirna_name_to_idx.values(al.target));
		for p = pre_idx'
			pre_expr.mean(p, s) = pre_expr.mean(p, s) + 1;
		end
	end
end

