
function expr = exon_expression_rnaseq(alignments, annotations)

global organism;

tmp = temporary('exon_expression_rnaseq');
bedtools_out = [tmp 'exon_expr.count'];

S = length(alignments.url);

% Count number of columns in annotations file.
data = readtable(annotations, 'Header', false, 'NumLines', 1);
num_annot_cols = length(data);

expr = struct;
expr.meta = alignments.meta;
expr.meta.type = 'Exon expression';

bam_urls = '';
for s = 1:S, bam_urls = [bam_urls alignments.url{s} '.bam ']; end

[status, out] = unix(sprintf( ...
	'bedtools multicov -bams %s -bed %s > %s', ...
	bam_urls, annotations, bedtools_out));
if status ~= 0, error('Bedtools coverage failed:\n%s', out); end

data = readtable(bedtools_out, 'Header', false, ...
	'Numeric', num_annot_cols+1:num_annot_cols+S);
if rx(annotations, '.*\.gtf')
	expr.rows.transcript = data{9};
	expr.rows.chromosome = chromosome_sym2num(data{1});
	expr.rows.position = [str2double(data{4}), str2double(data{5})];
elseif rx(annotations, '.*\.bed')
	expr.rows.gene = data{4};
	expr.rows.chromosome = chromosome_sym2num(data{1});
	expr.rows.position = [str2double(data{2}), str2double(data{3})];
end

expr.mean = nan(length(data{1}), S);
for s = 1:S, expr.mean(:, s) = data{num_annot_cols+s}; end

