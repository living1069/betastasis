classdef Lazy < handle
	properties (Hidden = true, Access = protected)
		Filename
		Done = false
	end
	
	methods
		function ret = lazy_get(obj, field)
			global organism;
			if ~obj.Done
				g = load([ppath '/organisms/' flatten_str(organism.Name) ...
					'/' flatten_str(organism.Version) ...
					'/' flatten_str(obj.Filename)]);
				eval(['organism.' obj.Filename ' = g.' ...
					flatten_str(obj.Filename) ';']);
			end
			obj.Done = true;
			eval(['ret = organism.' obj.Filename '.' field ';']);
		end

		%function obj = Lazy(file, fields)
		%	obj.Filename = file;
		%	for k = 1:length(fields)
		%		dp = addprop(obj, fields{k});
		%		dp.GetMethod = @lazy_get;
		%		dp.GetObservable = true;
		%		addlistener(obj, fields{k}, 'PreGet', @Lazy.get_listener);
		%	end
		%end
	end
end






