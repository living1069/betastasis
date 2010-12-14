classdef DiseaseOntology
	properties (GetAccess = public, SetAccess = private)
		Version
		ID
		AltID
		Name
		Description
		Obsolete
		IsA
		Synonym
		Subset
	end
	
	methods
	
	function idx = term_idx(obj, q)
		if ischar(q)
			tokens = regexpi(q, '^DOID:\s*(\d+)$', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1};
				id = str2double(tokens{1});
				idx = find(strcmp(id, obj.ID));
				if length(idx) == 1, return, end
			end
			
			id = str2double(q);
			if ~isnan(id)
				idx = find(strcmp(id, obj.ID));
				if length(idx) == 1, return, end
			end
			
			idx = find(strcmpi(q, obj.Name));
			if length(idx) == 1, return, end
		elseif isnumeric(q)
			idx = find(strcmp(q, obj.ID));
			if length(idx) == 1, return, end
		end
		error 'Invalid term identifier.';
	end
	
	function children = descendants(obj, term)
		children = [obj.term_idx(term)];
		k = 1;
		while k <= length(children)
			idx = children(k);
			children = [children obj.IsA(obj.IsA(:, 2) == idx, 1)'];
			k = k + 1;
		end
	
		children = unique(children(2:end));
	end

	function ontologies = DiseaseOntology(filepath, version)
	
		ontologies.Version = version;

		fid = fopen(filepath);

		prealloc = 100000;

		ontologies.ID = zeros(prealloc, 1);
		ontologies.AltID = cell(prealloc, 1);
		ontologies.Name = cell(prealloc, 1);
		ontologies.Obsolete = false(prealloc, 1);
		ontologies.Description = cell(prealloc, 1);
		ontologies.IsA = zeros(prealloc, 2);
		ontologies.Synonym = cell(prealloc, 1);
		ontologies.Subset = cell(prealloc, 1);

		parsing_term = false;
		id = 0;
		isa_count = 0;
		
		idx = 0;
		
		id_to_idx = containers.Map('KeyType', 'double', 'ValueType', 'double');

		while 1
			line = fgetl(fid);
			if line == -1, break, end
				
			if isempty(line), continue, end
			
			if line(1) == '['
				if strcmp(line, '[Term]')
					parsing_term = true;
				else
					parsing_term = false;
				end
				continue;
			end
			
			if parsing_term == false, continue, end
			
			tokens = regexp(line, '^id: DOID:(\d+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; id = str2double(tokens{1});
				idx = idx + 1;
				ontologies.ID(idx) = id;
				id_to_idx(id) = idx;
				continue;
			end
			
			tokens = regexp(line, 'alt_id: DOID:(\d+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; alt_id = str2double(tokens{1});
				ontologies.AltID{idx}(end+1) = alt_id;
			end
			
			tokens = regexp(line, 'name: (.+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; name = tokens{1};
				ontologies.Name{idx} = name;
				continue;
			end
			
			tokens = regexp(line, 'def: "(.+?)"', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; desc = tokens{1};
				ontologies.Description{idx} = desc;
				continue;
			end
			
			tokens = regexp(line, 'is_a: DOID:(\d+).*', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; base_term = str2double(tokens{1});
				isa_count = isa_count + 1;
				ontologies.IsA(isa_count, :) = [id base_term];
				continue;
			end
			
			tokens = regexp(line, 'synonym: (.+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; synonym = tokens{1};
				ontologies.Synonym{idx}{end+1} = synonym;
				continue;
			end

			tokens = regexp(line, 'subset: (.+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; subset = tokens{1};
				ontologies.Subset{idx}{end+1} = subset;
				continue;
			end
			
			if strcmp(line, 'is_obsolete: true')
				ontologies.Obsolete(idx) = true;
				continue;
			end
		end
		
		fclose(fid);
		
		ontologies.ID = ontologies.ID(1:idx);
		ontologies.AltID = ontologies.AltID(1:idx);
		ontologies.Name = ontologies.Name(1:idx);
		ontologies.Obsolete = ontologies.Obsolete(1:idx);
		ontologies.Description = ontologies.Description(1:idx);
		ontologies.IsA = ontologies.IsA(1:isa_count, :);
		ontologies.Synonym = ...
			ontologies.Synonym(1:idx);
		ontologies.Subset = ...
			ontologies.Subset(1:idx);
		ontologies.IsA(:) = cell2mat( ...
			id_to_idx.values(num2cell(ontologies.IsA(:))));
	end
	
	end
end


