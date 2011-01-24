
% ARACNE     Infer regulatory links between genes using the ARACNE algorithm
%
%    LINKS = ARACNE(EXPR) infers regulatory links between the genes available
%    in the microarray gene expression data EXPR. By default, links are inferred
%    between every pair of genes, with significance level 0.05 and
%    DPI tolerance 1.0.
%    
%    ARACNE(..., 'IncludeGenes', GENES) filters the ARACNE input expression
%    data so that only the genes specified by the index vector GENES are
%    included. This is not the same as subnetwork inference, where the input
%    data is not filtered at all. The indices in GENES must refer to the
%    currently selected organism.
%
%    ARACNE(..., 'Subnetwork', GENES) tells ARACNE to only infer regulatory
%    links for the genes named in cell array GENES. This means that
%    mutual information is not calculated for gene pairs where both genes
%    come from outside the subnetwork.
%
%    ARACNE(..., 'TFs', GENES) tells ARACNE to assume that only the genes
%    specified by the index vector GENES can act as transcription factors.
%    This affects the DPI pruning step.
%
%    ARACNE(..., 'ToleranceDPI', TOL) customizes the tolerance level TOL used
%    during the DPI pruning step. Default is 1.0 (no DPI pruning). Margolin
%    et al recommend values between 0.0 and 0.2 if DPI pruning is desired.
%
%    ARACNE(..., 'MaxThreads', THREADS) specifies the maximum number of
%    parallel ARACNE runs that this function is allowed to spawn. Defaults
%    to the maximum number of threads specified in the pipeline configuration.
%
%    ARACNE(..., 'Significance', ALPHA) specifies the significance level ALPHA
%    that ARACNE should use when pruning regulatory links between genes. No
%    familywise error rate correction is applied to the requested significance
%    level, so these must be calculated before calling this function.
%    Default is 1.0 (no significance testing, all MI values reported).

% Author: Matti Annala <matti.annala@tut.fi>

function links = aracne(expr, varargin)

global pipeline_config;
global organism;

max_threads = pipeline_config.MaxThreads;
include_genes = [];
subnetwork = [];
significance = 1.0;
dpi_tolerance = 1.0;
tfs = [];

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'IncludeGenes')
		include_genes = varargin{k+1};
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
	
	if strcmpi(varargin{k}, 'MaxThreads')
		max_threads = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Significance')
		significance = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

if isempty(subnetwork)
	subnetwork = input_genes;
else
	if iscellstr(subnetwork)
		found = ismember(subnetwork, organism.Genes.Name);
		not_found_idx = find(~found);
		for k = 1:length(not_found_idx)
			fprintf(1, 'Subnetwork gene %s is not recognized.\n', ...
				organism.Genes.Name{not_found_idx(k)});
		end
		sub_names = subnetwork(found);
		
		subnetwork = false(length(organism.Genes.Name), 1);
		[~, loc] = ismember(sub_names, organism.Genes.Name);
		subnetwork(loc) = true;
	elseif islogical(subnetwork)
		if length(subnetwork) ~= length(organism.Genes.Name)
			error(['The logical vector specifying the subnetwork must ' ...
				'have the same size as organism.Genes.Name.']);
		end
	else
		error(['Subnetwork must be specified either as a cell array of gene '...
			'names, or as a logical vector.']);
	end
end

input_genes = ~any(isnan(expr.Mean), 2);
if ~isempty(include_genes)
	input_genes = (input_genes & include_genes);
end

% Drop and warn about genes for which no expression values are available.
if ~all(subnetwork)
	no_expr = find(subnetwork & ~input_genes);
	for k = 1:length(no_expr)
		fprintf(1, 'No expression values for gene %s.\n', ...
			organism.Genes.Name{no_expr(k)});
	end
end

subnetwork = find(subnetwork & input_genes);


max_threads = min(max_threads, length(subnetwork));
parsize = floor(length(subnetwork) / max_threads);



fprintf(1, 'Constructing ARACNE gene expression input file...\n');
tmp = ptemp;
input_file = [tmp '.input_expr'];
fid = fopen(input_file, 'W');

fprintf(fid, 'Gene\tAnnotation');
for k = 1:length(expr.Meta.Sample.ID)
	fprintf(fid, '\t%s', expr.Meta.Sample.ID{k});
end
fprintf(fid, '\n');

input_genes = find(input_genes);
for g = 1:length(input_genes)
	gene = input_genes(g);
	fprintf(fid, '%s\t---', organism.Genes.Name{gene});
	fprintf(fid, '\t%f', expr.Mean(gene, :));
	fprintf(fid, '\n');
end
fclose(fid);

if max_threads > 1
	fprintf(1, 'Calculating mutual information (%d threads)...\n',max_threads);

	parnet_files = {};
	adj_files = {};

	pos = 1;
	for t = 1:max_threads
		if t < max_threads
			parnet = subnetwork(pos:pos+parsize-1);
		elseif t == max_threads
			parnet = subnetwork(pos:end);
		end
		
		parnet_files{t} = sprintf('%s.%d.parnet', tmp, t);
		fid = fopen(parnet_files{t}, 'W');
		for g = 1:length(parnet)
			fprintf(fid, '%s\n', organism.Genes.Name{parnet(g)});
		end
		fclose(fid);

		adj_files{t} = sprintf('%s.%d.adj', tmp, t);

		status = unix(sprintf(['%s/tools/aracne/aracne2 -H %s/tools/aracne ' ...
			'-s %s -i %s -o %s && touch %s.done &'], ppath, ppath, ...
			parnet_files{t}, input_file, adj_files{t}, adj_files{t}));
		if status ~= 0, error 'ARACNE network inference failed.'; end
		
		pos = pos + parsize;
	end

	for t = 1:max_threads
		while exist([adj_files{t} '.done']) == 0
			pause(1);
		end
		delete(parnet_files{t});
	end

	fprintf(1, 'Inferring network based on mutual information...\n');

	combined_file = [tmp '.total'];
	out_fid = fopen(combined_file, 'W');
	for t = 1:max_threads
		fid = fopen(adj_files{t});
		while 1
			line = fgetl(fid);
			if line == -1, break, end
			if t == 1 || line(1) ~= '>'
				fprintf(out_fid, '%s\n', line);
			end
		end
		fclose(fid);
	end
	fclose(out_fid);

	for t = 1:max_threads
		delete(adj_files{t});
		delete([adj_files{t} '.done']);
	end


	flags = '';

	if ~isempty(tfs)
		if islogical(tfs), tfs = find(tfs); end
		
		if iscellstr(tfs)
			found = ismember(tfs, organism.Genes.Name);
			not_found_idx = find(~found);
			for k = 1:length(not_found_idx)
				fprintf(1, 'Transcription factor %s is not recognized.\n', ...
					organism.Genes.Name{not_found_idx(k)});
			end
			tfs = tfs(found);
			
			[~, tfs] = ismember(tfs, organism.Genes.Name);
		end

		tf_file = [tmp '.tfs'];
		fid = fopen(tf_file, 'W');
		for g = 1:length(tfs)
			fprintf(fid, '%s\n', organism.Genes.Name{tfs(g)});
		end
		fclose(fid);
		flags = sprintf('%s -l %s', flags, tf_file);
	end


	final_file = [tmp '.final'];
	subnet_file = [tmp '.subnet'];

	fid = fopen(subnet_file, 'W');
	for g = 1:length(subnetwork)
		fprintf(fid, '%s\n', organism.Genes.Name{subnetwork(g)});
	end
	fclose(fid);



	status = unix(sprintf(['%s/tools/aracne/aracne2 -H %s/tools/aracne ' ...
		'%s -p %f -e %f -s %s -i %s -j %s -o %s'], ppath, ppath, flags, ...
		significance, dpi_tolerance, ...
		subnet_file, input_file, combined_file, final_file));
	if status ~= 0, error 'Final network construction failed.'; end


	
else
	fprintf(1, 'Inferring network (non-threaded)...\n');
	
	flags = '';

	if ~isempty(tfs)
		if islogical(tfs), tfs = find(tfs); end
		
		if iscellstr(tfs)
			found = ismember(tfs, organism.Genes.Name);
			not_found_idx = find(~found);
			for k = 1:length(not_found_idx)
				fprintf(1, 'Transcription factor %s is not recognized.\n', ...
					organism.Genes.Name{not_found_idx(k)});
			end
			tfs = tfs(found);
			
			[~, tfs] = ismember(tfs, organism.Genes.Name);
		end

		tf_file = [tmp '.tfs'];
		fid = fopen(tf_file, 'W');
		for g = 1:length(tfs)
			fprintf(fid, '%s\n', organism.Genes.Name{tfs(g)});
		end
		fclose(fid);
		flags = sprintf('%s -l %s', flags, tf_file);
	end

	final_file = [tmp '.final'];
	subnet_file = [tmp '.subnet'];

	fid = fopen(subnet_file, 'W');
	for g = 1:length(subnetwork)
		fprintf(fid, '%s\n', organism.Genes.Name{subnetwork(g)});
	end
	fclose(fid);

	status = unix(sprintf(['%s/tools/aracne/aracne2 -H %s/tools/aracne ' ...
		'%s -p %f -e %f -s %s -i %s -o %s'], ppath, ppath, flags, ...
		significance, dpi_tolerance, subnet_file, input_file, final_file));
	if status ~= 0, error 'Final network construction failed.'; end
end

links = struct;
links.Genes = [];
links.MI = [];

fid = fopen(final_file);
while 1
	line = fgetl(fid);
	if line == -1, break, end
		
	data = textscan(line, '%s', -1, 'Delimiter', '\t');
	data = data{1};

	regulator = gene_idx(data{1});
	for k = 2:2:length(data)
		links.Genes(end+1, :) = [regulator gene_idx(data{k})];
		links.MI(end+1, 1) = str2double(data{k+1});
	end
end
fclose(fid);

%if nargout == 0
%	for k = 1:length(links.MI)
%		fprintf(1, '%s - %s: %.4f\n', organism.Genes.Name{links.Genes(k,1)}, ...
%			organism.Genes.Name{links.Genes(k, 2)}, links.MI(k));
%	end
%end

