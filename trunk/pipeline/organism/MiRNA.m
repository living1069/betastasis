classdef MiRNA < Lazy
	properties
		Name
		Accession
		Sequence
	end
	
	methods
		function obj = MiRNA(), obj.Filename = 'miRNA'; end
		
		function r = get.Name(obj), r = obj.lazy_get('Name'); end
		function r = get.Accession(obj), r = obj.lazy_get('Accession'); end
		function r = get.Sequence(obj), r = obj.lazy_get('Sequence'); end
	end
end






