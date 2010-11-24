classdef HttpsRepository < handle
	properties (Access = public)
		Name = ''
		URL = ''
		Username = ''
		CheckCertificate = false
	end
	
	properties (GetAccess = private, SetAccess = public, Hidden = true)
		Password = ''
	end
	
	properties (Access = private, Hidden = true)
		SampleCache = ''
	end

	methods
		function obj = HttpsRepository(name, url, user)
			obj.Name = name;
			obj.URL = url;
			obj.Username = user;
			obj.SampleCache = ptemp;
		end
		
		function items = contents(obj)
			items = {};
			
			obj.prompt_password();
			
			if isunix || ismac
				out_dir = ptemp;
				obj.download([obj.URL '/index.txt'], out_dir);
			else
				error 'HTTPS repositories are not supported on Windows.';
			end
			
			fid = fopen([out_dir '/index.txt']);
			paths = textscan(fid, '%s', 'Delimiter', '\t');
			paths = paths{1};
			fclose(fid);
			
			rmdir(out_dir, 's');
			
			for k = 1:length(paths)
				slashes = find(paths{k} == '/');
				for s = 1:length(slashes)
					items{end + 1, 1} = paths{k}(1:slashes(s));
				end
				items{end + 1, 1} = paths{k};
			end
			
			items = unique(items);
		end

		function [] = cache(obj, files)
			dirs = regexprep(files, '(.*)/.*', '$1');
			[ds, ~, ds_labels] = unique(dirs);
			
			obj.prompt_password();
			fprintf(1, '\b\b\b(downloading)');
			
			if isunix || ismac
				for k = 1:length(ds)
					out_dir = [obj.SampleCache '/' ds{k}];
					ds_urls = files(ds_labels == k);
					obj.download(strcat(obj.URL, '/', ds_urls), out_dir);
				end
			else
				error 'HTTPS repositories are not supported on Windows.';
			end
			
			fprintf(1, [repmat('\b', 1, 13) repmat(' ', 1, 13) ...
				repmat('\b', 1, 13) '00%%']);
		end
		
		function data = load_cached(obj, file)
			global pipeline_config;
			cache_path = [obj.SampleCache '/' file];
			data = load(cache_path);
			delete(cache_path);
		end
		
		function data = load(obj, path)
			obj.prompt_password();
			
			if isunix || ismac
				out_dir = ptemp;
				obj.download([obj.URL '/' path], out_dir);
			else
				error 'HTTPS repositories are not supported on Windows.';
			end
			
			fid = fopen([out_dir '/index.txt']);
			data = load([out_dir '/' path_strip_dir(path)]);
		end
		
		function [] = prompt_password(obj)
			if isempty(obj.Password) || ~ischar(obj.Password)
				if ispc
					obj.Password = passwordEntryDialog('WindowName', ...
						sprintf('Password for user %s at %s:', obj.Username, ...
						obj.URL), 'CheckPasswordLength', false);
					if obj.Password == -1, obj.Password = ''; end
					if isempty(obj.Password)
						error(['You must specify a password in ' ...
							'order to use an HTTPS repository.']);
					end
				else				
					obj.Password = input('Password: ', 's');
					ccount = length('Password: ') + length(obj.Password)+1;
					fprintf(1, [repmat('\b', 1, ccount) ...
						repmat(' ', 1, ccount) '\n']);
				end
			end
		end
		
		function [] = download(obj, urls, out_dir)
			if ischar(urls), urls = { urls }; end
				
			flags = '';
			if ~obj.CheckCertificate
				flags = [flags ' --no-check-certificate'];
			end
			
			cmd_file = ptemp;
			fid = fopen(cmd_file, 'W');
			for k = 1:length(urls)
				url = regexprep(urls{k}, '(https?://)(.*)$', ...
					['$1' obj.Username ':' obj.Password '@$2']);
				fprintf(fid, '%s\n', url);
			end
			fclose(fid);
			
			[~, ~] = mkdir(out_dir);
			[status, output] = unix(['wget ' flags ' -i ' cmd_file ' -P ' ...
				out_dir]);
			delete(cmd_file);
			
			if status ~= 0
				error('Download failed:\n%s\n', output);
			end
		end
	end
end

