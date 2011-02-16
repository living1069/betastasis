
% ARACNE     Infer regulatory links between genes using the ARACNE algorithm
%
%    LINKS = ARACNE(EXPR) infers regulatory links between the genes available
%    in the microarray gene expression data EXPR. By default, links are inferred
%    between every pair of genes, with significance level 0.05 and
%    DPI tolerance 1.0.
%
%    The function can also be used to infer regulatory links between small RNA
%    and genes. This is done by concatenating the expression matrices into one
%    large matrix, which is then given as the first argument to ARACNE(). The
%    names of the genetic features must then be provided using the optional
%    argument 'FeatureNames' (see below).
%
%    Any expression values provided to this function should be given in the
%    natural scale: a log-2 transformation will be applied to all expression
%    values before the ARACNE executable is run.
%
%    All expression features with even a single NaN value in one of the samples
%    are discarded and not used in the analysis. If you wish to intentionally
%    use only a subset of all available expression values, set the rest to a
%    NaN value.
%
%    ARACNE(..., 'Subnetwork', GENES) tells ARACNE to only infer regulatory
%    links for the genes identified by the logical vector ROWS. This means that
%    mutual information is not calculated for gene pairs where both genes are
%    outside the subnetwork.
%
%    ARACNE(..., 'TFs', GENES) tells ARACNE to assume that only the genes
%    specified by the index vector GENES can act as transcription factors.
%    This affects the DPI pruning step.
%
%    ARACNE(..., 'ToleranceDPI', TOL) customizes the tolerance level TOL used
%    during the DPI pruning step. Default is 1.0 (no DPI pruning). Margolin
%    et al recommend values between 0.0 and 0.2 if DPI pruning is desired.
%
%    ARACNE(..., 'Significance', ALPHA) specifies the significance level ALPHA
%    that ARACNE should use when pruning regulatory links between features. No
%    familywise error rate correction is applied to the requested significance
%    level, so these must be calculated before calling this function.
%    Default is 1.0 (no significance testing, all MI values reported).

% Author: Matti Annala <matti.annala@tut.fi>

function [links, feature_names] = aracne(expr, varargin)

global organism;
genes = organism.Genes;

max_threads = 1;
feature_names = [];
subnetwork = [];
significance = 1.0;
dpi_tolerance = 1.0;
tfs = [];

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'FeatureNames')
		feature_names = varargin{k+1};
		if ~iscellstr(feature_names)
			error 'Feature names must be provided as a cell array of strings.';
		end
		continue;
	end
	
	if strcmpi(varargin{k}, 'Subnetwork')
		subnetwork = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'TFs')
		tfs = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'ToleranceDPI')
		dpi_tolerance = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Significance')
		significance = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if isnumeric(expr)
	log_expr = log2(expr);
elseif isfield(expr, 'Mean')
	log_expr = log2(expr.Mean),
else
	error 'Expression data was provided in an unrecognized format.';
end

G = size(log_expr, 1);
input_genes = ~any(isnan(log_expr), 2);

if isempty(subnetwork)
	subnetwork = input_genes;
else
	if islogical(subnetwork)
		if length(subnetwork) ~= G
			error(['The subnetwork logical vector must have as many ' ...
				'elements as there are rows in the expression matrix.']);
		end
	else
		error 'Subnetwork must be specified as a logical vector.';
	end
end

if ~isempty(feature_names)
	if length(feature_names) ~= G
		error(['The number of feature names provided must equal the number ' ...
			'of rows in the expression matrix.']);
	end
else
	if G ~= length(genes.Name)
		error(['No feature names were provided, and the expression matrix ' ...
			'has more rows than there are genes for the organism.']);
	end
	feature_names = genes.Name;
end
	
	
% Drop and warn about genes for which no expression values are available.
if ~all(subnetwork)
	no_expr = find(subnetwork & ~input_genes);
	for k = 1:length(no_expr)
		fprintf(1, 'No expression values for %s.\n', ...
			feature_names{no_expr(k)});
	end
end

subnetwork = find(subnetwork & input_genes);


%max_threads = min(max_threads, length(subnetwork));
%parsize = floor(length(subnetwork) / max_threads);



fprintf(1, 'Constructing ARACNE gene expression input file...\n');
tmp = ptemp;
input_file = [tmp '.input_expr'];
fid = fopen(input_file, 'W');

fprintf(fid, 'Gene\tAnnotation');
for k = 1:size(log_expr, 2)
	fprintf(fid, '\tSample_%d', k);
end
fprintf(fid, '\n');

input_genes = find(input_genes);
for g = 1:length(input_genes)
	gene = input_genes(g);
	fprintf(fid, '%s\t---', feature_names{gene});
	fprintf(fid, '\t%f', log_expr(gene, :));
	fprintf(fid, '\n');
end
fclose(fid);



fprintf(1, 'Inferring network...\n');

flags = '';

if ~isempty(tfs)
	if islogical(tfs), tfs = find(tfs); end
	
	if iscellstr(tfs)
		found = ismember(tfs, feature_names);
		not_found_idx = find(~found);
		for k = 1:length(not_found_idx)
			fprintf(1, 'Transcription factor %s is not recognized.\n', ...
				feature_names{not_found_idx(k)});
		end
		tfs = tfs(found);
		
		[~, tfs] = ismember(tfs, feature_names);
	end

	tf_file = [tmp '.tfs'];
	fid = fopen(tf_file, 'W');
	for g = 1:length(tfs)
		fprintf(fid, '%s\n', feature_names{tfs(g)});
	end
	fclose(fid);
	flags = sprintf('%s -l %s', flags, tf_file);
end

final_file = [tmp '.final'];
subnet_file = [tmp '.subnet'];

fid = fopen(subnet_file, 'W');
for g = 1:length(subnetwork)
	fprintf(fid, '%s\n', feature_names{subnetwork(g)});
end
fclose(fid);

status = unix(sprintf(['%s/tools/aracne/aracne2 -H %s/tools/aracne ' ...
	'%s -p %f -e %f -s %s -i %s -o %s'], ppath, ppath, flags, ...
	significance, dpi_tolerance, subnet_file, input_file, final_file));
if status ~= 0, error 'Final network construction failed.'; end

	
	
links = struct;
links.Genes = zeros(1e6, 2);
links.MI = zeros(1e6, 1);
L = 0;

fid = fopen(final_file);
while 1
	line = fgetl(fid);
	if line == -1, break, end
		
	data = textscan(line, '%s', -1, 'Delimiter', '\t');
	data = data{1};

	regulator = find(strcmp(data{k}, feature_names));
	for k = 2:2:length(data)
		L = L + 1;
		links.Genes(L, :) = [regulator find(strcmp(data{k}, feature_names))];
		links.MI(L, 1) = str2double(data{k+1});
	end
end
fclose(fid);

links.Genes = links.Genes(1:L, :);
links.MI = links.MI(1:L);

