classdef FilePool < handle
	properties (GetAccess = public, SetAccess = private)
		Paths = {}
	end
	
	properties (Access = private)
		Temporary = []
		TempPrefix = ''
	end
	
	methods
		function obj = FilePool(paths)
			if nargin == 1
				obj.Paths = paths;
				obj.Temporary = false(1, length(files));
			end
			obj.TempPrefix = ptemp;
		end
		
		function path = static(obj, path)
			obj.Paths{end+1, 1} = path;
			obj.Temporary(1, end+1) = false;
		end
		
		function path = temp(obj, suffix)
			path = [obj.TempPrefix '.' suffix];
			obj.Paths{end+1, 1} = path;
			obj.Temporary(1, end+1) = true;
		end
		
		function len = length(obj)
			len = length(obj.Paths);
		end
		
		function ret = subsref(obj, s)
			if s(1).type == '.'
				ret = builtin('subsref', obj, s);
				return;
			end
		
			if length(s) ~= 1, error 'Invalid access.'; end
			ret = subsref(obj.Paths, s);
		end
		
		function delete(obj)
			for f = 1:length(obj.Paths)
				if obj.Temporary(f), safe_delete(obj.Paths{f}); end
			end
		end
	end
end
	

