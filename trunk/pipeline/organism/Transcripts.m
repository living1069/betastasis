classdef Transcripts < Lazy
	properties
		Name
		Gene
		Sequence
		CDS
	end
	
	methods
		function obj = Transcripts(), obj.Filename = 'Transcripts'; end
		
		function r = get.Name(obj), r = obj.lazy_get('Name'); end
		function r = get.Gene(obj), r = obj.lazy_get('Gene'); end
		function r = get.Sequence(obj), r = obj.lazy_get('Sequence'); end
		function r = get.CDS(obj), r = obj.lazy_get('CDS'); end
	end
end






