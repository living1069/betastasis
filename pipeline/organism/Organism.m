classdef Organism < dynamicprops
	properties (SetAccess = private)
		Name
		Version
		Chromosomes
		Genes
		Transcripts
		Exons
		miRNA
		pre_miRNA
		SNPs
		Cells
		Ontologies
	end
	
	properties (Hidden = true, Access = private)
		Fields
	end
	
	methods
		function obj = Organism(name, version)
			obj.Name = name;
			obj.Name(1) = upper(name(1));
			obj.Version = version;
			obj.Fields = containers.Map;
			
			obj.Fields('Name') = obj.Name;
			obj.Fields('Version') = obj.Version;
			
			obj.add_field('Chromosomes');
			obj.add_field('Genes');
			obj.add_field('Transcripts');
			obj.add_field('Exons');
			obj.add_field('miRNA');
			obj.add_field('pre_miRNA');
			obj.add_field('SNPs');
			obj.add_field('Cells');
			obj.add_field('Ontologies.DARNED');
			obj.add_field('Ontologies.Disease');
			obj.add_field('Ontologies.GO');
			obj.add_field('Ontologies.OMIM');
			obj.add_field('Ontologies.PC');
		end
		
		function ret = add_field(obj, field)
			rel_path = flatten_str(strrep(field, '.', '/'));
			if ~exist([obj.org_path '/' rel_path '.mat']), return, end
				
			field = field;
			
			dots = [find(field == '.') length(field)+1];
			for k = 1:length(dots)
				base = field(1:dots(k)-1);
				if ~obj.Fields.isKey(base)
					if k == 1
						eval(['obj.' base ' = struct;']);
					else
						prev_base = field(1:dots(k-1)-1);
						prev_struct = obj.Fields(prev_base);
						eval(['prev_struct.' field(dots(k-1)+1:dots(k)-1) ...
							' = struct;']);
						obj.Fields(prev_base) = prev_struct;
					end
					
					if k == length(dots)
						obj.Fields(base) = [];
					else
						obj.Fields(base) = struct;
					end
				end
			end
		end
	
		function ret = subsref(obj, s)
			last_match = '';
			field = '';
			for k = 1:length(s)+1
				if k == length(s)+1, break, end
				
				if ~strcmp(s(k).type, '.')
					if isempty(last_match), error 'Invalid access.'; end
					break;
				end
				
				field = [field '.' s(k).subs];
				if ~obj.Fields.isKey(field(2:end))
					if isempty(last_match), error 'Invalid access.'; end
					break;
				end
				
				last_match = field(2:end);
			end
			
			field = last_match;
			
			ret = obj.Fields(field);
			if isempty(ret)
				rel_path = flatten_str(strrep(field, '.', '/'));
				ret = load([obj.org_path '/' rel_path]);
				f = fieldnames(ret);
				eval(['ret = ret.' f{1} ';']);
				obj.Fields(field) = ret;
			end
			
			if k ~= length(s)+1
				ret = subsref(ret, s(k:end));
			end
		end
		
		function ret = org_path(obj)
			ret = [ppath '/organisms/' flatten_str(obj.Name) ...
				'/' flatten_str(obj.Version)];
		end
	end
end





