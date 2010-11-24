classdef deferred < handle
	properties (Hidden = true, Access = private)
		DataPath = '';
		Listeners = {};
	end

	methods
		function obj = deferred(path)
			prop = properties(obj);
			for k = 1:length(prop)
				obj.Listeners{end + 1} = addlistener(obj, prop{k}, ...
					'PreGet', @obj.lazy_preinit);
			end
			obj.DataPath = path;
		end
		
		function lazy_preinit(obj, src, event)
			data = load(obj.DataPath);
			fields = fieldnames(data);
						
			eval(['data = data.' fields{1} ';']);
			fields = fieldnames(data);
			for k = 1:length(fields)
				eval(['obj.' fields{k} ' = data.' fields{k} ';']);
			end
			for k = 1:length(obj.Listeners), delete(obj.Listeners{k}); end
		end
	end
end

