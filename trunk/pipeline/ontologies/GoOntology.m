classdef GoOntology
	properties (GetAccess = public, SetAccess = private)
		Name
		AltID
		Namespace
		Description
		Obsolete
		IsA
		Regulates
		PartOf
	end
	
	properties (Hidden = true)
		TermGenes
	end
	
	methods
	
	function idx = term_idx(obj, id)
		if ischar(id)
			tokens = regexpi(id, '^GO:\s*(\d+)$', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1};
				idx = str2double(tokens{1});
				return;
			end
			
			idx = str2double(id);
			if ~isnan(idx), return, end
			
			idx = find(strcmpi(id, obj.Name));
			if length(idx) == 1, return, end
		elseif isnumeric(id)
			idx = id;
			return;
		end
		error 'Invalid GO term identifier.';
	end
	
	function genes = genes(obj, term)
		genes = [];
		term_stack = [obj.term_idx(term)];
		while ~isempty(term_stack)
			idx = term_stack(1);
			term_stack = term_stack(2:end);
			genes = [genes obj.TermGenes{idx}];
			
			term_stack = [term_stack obj.IsA(obj.IsA(:, 2) == idx, 1)'];
			term_stack = [term_stack obj.PartOf(obj.PartOf(:, 2) == idx, 1)'];
			term_stack = unique(term_stack);
		end
		
		genes = unique(genes);
	end

	function ontologies = GoOntology(ontologies)

		ontology_file = '';
		onto_assoc_file = '';

		files = dir('.');
		for k = 1:length(files)
			if files(k).isdir, continue, end
				
			if regexpi(files(k).name, 'gene_ontology\..+\.obo')
				if ~isempty(ontology_file)
					error 'Found more than one ontology description file.';
				end
				ontology_file = files(k).name;
				fprintf(1, 'Reading ontology information from %s...\n', ...
					ontology_file);
			end
			
			if regexpi(files(k).name, 'gene_association\.goa_human')
				onto_assoc_file = files(k).name;
				fprintf(1, 'Reading gene associations from %s...\n', ...
					onto_assoc_file);
			end
		end

		if isempty(ontology_file)
			error 'No ontology description file found.';
		end

		if isempty(onto_assoc_file)
			error 'No gene ontology association file found.';
		end

		fid = fopen(ontology_file);

		prealloc = 100000;

		ontologies.AltID = cell(prealloc, 1);
		ontologies.Name = cell(prealloc, 1);
		ontologies.Namespace = repmat('-', prealloc, 1);
		ontologies.Obsolete = false(prealloc, 1);
		ontologies.Description = cell(prealloc, 1);
		ontologies.IsA = zeros(prealloc, 2);
		ontologies.Regulates = zeros(prealloc, 2);
		ontologies.PartOf = zeros(prealloc, 2);

		parsing_term = false;
		id = 0;
		highest_id = 0;
		isa_count = 0;
		regulates_count = 0;
		partof_count = 0;

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
			
			tokens = regexp(line, '^id: GO:(\d+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; id = str2double(tokens{1});
				highest_id = max(id, highest_id);
				continue;
			end
			
			tokens = regexp(line, 'alt_id: GO:(\d+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; alt_id = str2double(tokens{1});
				ontologies.AltID{id}(end+1) = alt_id;
			end
			
			tokens = regexp(line, 'name: (.+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; name = tokens{1};
				ontologies.Name{id} = name;
				continue;
			end
			
			tokens = regexp(line, 'def: "(.+?)"', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; desc = tokens{1};
				ontologies.Description{id} = desc;
				continue;
			end
			
			tokens = regexp(line, 'namespace: (.+)$', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; ns = tokens{1};
				switch ns
				case 'cellular_component', ontologies.Namespace(id) = 'C';
				case 'biological_process', ontologies.Namespace(id) = 'P';
				case 'molecular_function', ontologies.Namespace(id) = 'F';
				end
				continue;
			end
			
			tokens = regexp(line, 'is_a: GO:(\d+).*', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; base_term = str2double(tokens{1});
				isa_count = isa_count + 1;
				ontologies.IsA(isa_count, :) = [id base_term];
				continue;
			end
			
			tokens = regexp(line, 'relationship: .*regulates GO:(\d+)', ...
				'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; regulatee = str2double(tokens{1});
				regulates_count = regulates_count + 1;
				ontologies.Regulates(regulates_count, :) = [id regulatee];
				continue;
			end
			
			tokens = regexp(line, 'relationship: part_of GO:(\d+)', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; whole = str2double(tokens{1});
				partof_count = partof_count + 1;
				ontologies.PartOf(partof_count, :) = [id whole];
				continue;
			end

			if strcmp(line, 'is_obsolete: true')
				ontologies.Obsolete(id) = true;
				continue;
			end
		end

		fclose(fid);

		ontologies.AltID = ontologies.AltID(1:highest_id);
		ontologies.Name = ontologies.Name(1:highest_id);
		ontologies.Namespace = ontologies.Namespace(1:highest_id);
		ontologies.Obsolete = ontologies.Obsolete(1:highest_id);
		ontologies.Description = ontologies.Description(1:highest_id);
		ontologies.IsA = ontologies.IsA(1:isa_count, :);
		ontologies.Regulates = ontologies.Regulates(1:regulates_count, :);
		ontologies.PartOf = ontologies.PartOf(1:partof_count, :);




		fid = fopen('gene_association.goa_human');

		ontologies.TermGenes = cell(highest_id, 1);

		while 1
			line = fgetl(fid);
			if line == -1, break, end
			
			tokens = regexp(line, 'UniProtKB\t\w+\t(\w+)\t\w*\tGO:(\d+).*', ...
				'tokens');
			if length(tokens) == 1
				tokens = tokens{1};
				gene = tokens{1};
				ontology = str2double(tokens{2});
				
				idx = gene_idx(gene);
				if ~isnan(idx)
					ontologies.TermGenes{ontology} = ...
						[ontologies.TermGenes{ontology} idx];
				end
			end
		end

		fclose(fid);
	end

	
	

	end
end


