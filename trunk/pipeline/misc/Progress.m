classdef Progress < handle
	properties (SetAccess = private)
		Percentage = 0;
	end
	
	methods
		function obj = Progress()
			fprintf(1, 'Progress: 00%%');
		end
		function update(obj, val)
			if floor(val * 100) > obj.Percentage
				obj.Percentage = floor(val * 100);
				fprintf(1, '\b\b\b%02d%%', obj.Percentage);
				if obj.Percentage == 100
					fprintf(1, '\n');
				end
			end
		end
	end
end

