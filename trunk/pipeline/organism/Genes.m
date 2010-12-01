classdef Genes < Lazy
	properties
		Name
		EntrezID
		TranscriptCount
		Transcripts
		Chromosome
		Strand
		Position
	end
	
	methods
		function obj = Genes(), obj.Filename = 'Genes'; end
		
		function r = get.Name(obj), r = obj.lazy_get('Name'); end
		function r = get.Chromosome(obj), r = obj.lazy_get('Chromosome'); end
		function r = get.Strand(obj), r = obj.lazy_get('Strand'); end
		function r = get.Position(obj), r = obj.lazy_get('Position'); end
		function r = get.EntrezID(obj), r = obj.lazy_get('EntrezID'); end
		function r = get.Transcripts(obj), r = obj.lazy_get('Transcripts'); end
		function r = get.TranscriptCount(obj)
			r = obj.lazy_get('TranscriptCount');
		end
	end
end






