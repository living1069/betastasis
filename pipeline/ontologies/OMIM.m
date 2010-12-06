classdef OMIM
	properties (GetAccess = public, SetAccess = private)
		Version
		Diseases
	end
	
	methods
	
	function obj = OMIM(filepath, version)
	
		global organism;
	
		obj.Version = version;
		obj.Diseases = cell(length(organism.Genes.Name), 1);

		fid = fopen(filepath);
		
		while 1
			line = fgetl(fid);
			if line == -1, break, end
			
			tokens = regexp(line, '\s*(?{.*?}|[^|])*?\s*(\||$)', 'tokens');
			for k = 1:length(tokens)
				token = tokens{k}; cols{k} = token{1};
			end
			
			gene_names = cols{6};
			tokens = regexp(gene_names, '([^\s,]+)\s*(,|$)', 'tokens');
			if length(tokens) == 0, error 'No gene names found.'; end
				
			gene_names = {};
			for k = 1:length(tokens)
				token = tokens{k}; gene_names{k} = token{1};
			end
			
			diseases = [cols{14} cols{15} cols{16}];
			
			tokens = regexp(diseases, '\s*(.+?)\s*(;|$)', 'tokens');
				
			diseases = {};
			for k = 1:length(tokens)
				token = tokens{k}; dis = token{1};
				dis = regexprep(dis, '{(.*)}', '$1');
				dis = regexprep(dis, '\[(.*)\]', '$1');
				dis = regexprep(dis, '\s*\(\d+\)\s*$', '');
				dis = regexprep(dis, '\s*(,\s*(\d+)\s*)+$', '');
				diseases{k} = dis;
			end
			
			if length(diseases) == 0, continue, end
			
			for k = 1:length(gene_names)
				idx = gene_idx(gene_names{k});
				if ~isnan(idx)
					obj.Diseases{idx} = diseases;
					break;
				end
				
				%if k == length(gene_names)
				%	fprintf(1, 'Cannot find gene.\n');
				%	gene_names
				%end
			end
		end
		
		fclose(fid);
	end
	
	end
end


