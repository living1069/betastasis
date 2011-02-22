classdef LocalRepository < handle
	properties (Access = public)
		Name
		URL
	end

	methods
		function obj = LocalRepository(name, url)
			obj.Name = name;
			obj.URL = url;
		end
	
		function items = contents(obj)
			items = {};
			recurse_stack = {''};

			% Recursively look for all dataset directories.
			while 1
				if isempty(recurse_stack), break, end
					
				root = recurse_stack{1};
				recurse_stack = recurse_stack(2:end, :);
				
				dirs = dir([obj.URL '/' root]);
				for k = 1:length(dirs)
					if dirs(k).name(1) == '.' || ~dirs(k).isdir, continue, end
					d = [root dirs(k).name];
					if exist([obj.URL '/' d '/metadata.mat'])
						items{end + 1, 1} = d;
					else
						items{end + 1, 1} = [d '/'];
						recurse_stack{end + 1, 1} = [d '/'];
					end
				end
			end

			items = sort(items);
		end
		
		function cached = cache(obj, files)
			cached = strcat([obj.URL '/'], files);
		end
		
		function [] = clear_cache(obj, cached)
			return;
		end
	end
end

