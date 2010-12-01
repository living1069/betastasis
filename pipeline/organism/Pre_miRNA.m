classdef Pre_miRNA < Lazy
	properties
		Name
		Accession
		Sequence
		MatureCount
		Matures
		MatureOffsets
	end
	
	methods
		function obj = Pre_miRNA(), obj.Filename = 'pre_miRNA'; end
		
		function r = get.Name(obj), r = obj.lazy_get('Name'); end
		function r = get.Accession(obj), r = obj.lazy_get('Accession'); end
		function r = get.Sequence(obj), r = obj.lazy_get('Sequence'); end
		function r = get.MatureCount(obj), r = obj.lazy_get('MatureCount'); end
		function r = get.Matures(obj), r = obj.lazy_get('Matures'); end
		function r = get.MatureOffsets(obj)
			r = obj.lazy_get('MatureOffsets');
		end
	end
end






