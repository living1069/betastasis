classdef Chromosomes < Lazy
	properties
		Name
		Length
	end
	
	methods
		function obj = Chromosomes(), obj.Filename = 'Chromosomes'; end
		
		function r = get.Name(obj), r = obj.lazy_get('Name'); end
		function r = get.Length(obj), r = obj.lazy_get('Length'); end
	end
end






