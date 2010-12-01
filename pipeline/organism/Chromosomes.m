classdef Chromosomes < Lazy
	properties
		Name
		Length
		Sequence
	end
	
	methods
		function obj = Chromosomes(), obj.Filename = 'Chromosomes'; end
		
		function r = get.Name(obj), r = obj.lazy_get('Name'); end
		function r = get.Length(obj), r = obj.lazy_get('Length'); end
		function r = get.Sequence(obj), r = obj.lazy_get('Sequence'); end
		
	end
end






