classdef PathwayCommons
	properties (GetAccess = public, SetAccess = private)
		Version
		Name
		Source
		Genes
	end
	
	methods
	
	function obj = PathwayCommons(filepath, version)
	
		global organism;
		
		entrez_to_idx = containers.Map(organism.Genes.EntrezID, ...
			num2cell(1:length(organism.Genes.EntrezID)));
	
		obj.Version = version;
		obj.Name = {};
		obj.Source = {};
		obj.Genes = {};

		fid = fopen(filepath);
		
		while 1
			line = fgetl(fid);
			if line == -1, break, end
			
			tabs = find(line == sprintf('\t'));
			obj.Name{end+1, 1} = line(1:tabs(1)-1);
			obj.Source{end+1, 1} = line(tabs(1)+1:tabs(2)-1);
			obj.Genes{end+1, 1} = [];
			
			tokens = regexp(line, '\t\d+:protein:.+?:.+?:.+?:(\d+)', 'tokens');
			if length(tokens) == 0, continue, end
			
			for k = 1:length(tokens)
				token = tokens{k}; entrez = str2double(token{1});
				if ~entrez_to_idx.isKey(entrez)
					fprintf(1, 'No match for entrez ID %d.\n', entrez);
					continue;
				end
				obj.Genes{end}(end+1) = entrez_to_idx(entrez);
			end
		end
		
		fclose(fid);
	end
	
	end
end


