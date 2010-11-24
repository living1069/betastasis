classdef SshRepository < handle
	properties (Access = public)
		Name = ''
		URL = ''
		Hostname = ''
		Path = ''
		Username = ''
	end
	
	properties (Access = private, Hidden = true)
		Password = ''
	end

	methods
		function obj = SshRepository(name, url)
			obj.Name = name;
			obj.URL = url;
			
			tokens = regexp(obj.URL, '^(.+?)@(.+?):(/.+)$', 'tokens');
			if length(tokens) ~= 1
				error(['Invalid SSH location specified. The location must ' ...
					'be specified in the format user@host:/absolute/path']);
			end
			
			tokens = tokens{1};
			obj.Username = tokens{1};
			obj.Hostname = tokens{2};
			obj.Path = tokens{3};
		end
		
		function items = contents(obj)
			items = {};
			
			obj.prompt_password();
			
			if isunix || ismac
				[status, output] = unix(['ssh -o BatchMode=true ' ...
					obj.Username '@' ...
					obj.Hostname ' find ' obj.Path ' -name metadata.mat']);
				if status ~= 0
					if regexp(output, 'Permission denied')
						fprintf(1,'Public key authentication failed.\n');
					else
						fprintf(1,'Could not retrieve dataset list via SSH:\n');
						fprintf(1,'%s\n', output);
					end
					return;
				end
			else
				[status, output] = dos(['"' ppath ...
					'/tools/windows/plink.exe" -pw ' obj.Password ' '...
					obj.Username '@' obj.Hostname ...
					' find ' obj.Path ' -name metadata.mat']);
				if status ~= 0
					fprintf(1,'Could not retrieve dataset list via SSH:\n');
					fprintf(1, '%s\n', output);
					return;
				end
			end
			
			paths = textscan(output, '%s', 'Delimiter', '\t');
			paths = paths{1};
			paths = strrep(paths, obj.Path, '');
			for k = 1:length(paths)
				if paths{k}(1) == '/', paths{k} = paths{k}(2:end); end
			end
			
			for k = 1:length(paths)
				slashes = find(paths{k} == '/');
				if length(slashes) == 0, continue, end
					
				for s = 1:length(slashes)-1
					items{end + 1, 1} = paths{k}(1:slashes(s));
				end
				items{end + 1, 1} = paths{k}(1:slashes(end)-1);
			end
			
			items = unique(items);
		end

		function [] = cache(obj, files)
			global pipeline_config;
			
			cache_tmp = ptemp;
			
			tar_cmd = ['cd ' obj.Path ' && tar -cf ' cache_tmp ' '];
			for k = 1:length(files)
				tar_cmd = [tar_cmd files{k} ' '];
			end
			
			fprintf(1, '\b\b\b(downloading)');
			
			obj.prompt_password();
			if isunix || ismac
				[status, output] = unix(sprintf( ...
					'ssh -o BatchMode=true %s@%s "%s"', ...
					obj.Username, obj.Hostname, tar_cmd));
				if status ~= 0
					fprintf(1,'\nCould not prepare dataset for transfer:\n');
					fprintf(1, '%s\n', output);
					return;
				end
				
				obj.prompt_password();
				[status, ~] = unix(sprintf( ...
					'scp -o BatchMode=true %s@%s:%s %s', ...
					obj.Username, obj.Hostname, cache_tmp, ...
					pipeline_config.TempDir));
				if status ~= 0, return, end
			else
				cmd_file = ptemp;
				fid = fopen(cmd_file, 'W');
				fprintf(fid, '%s\n', tar_cmd);
				fclose(fid);
				
				[status, output] = dos(['"' ppath ...
					'/tools/windows/plink.exe" -m ' cmd_file ...
					' -pw ' obj.Password ' ' obj.Username '@' obj.Hostname]);
				if status ~= 0
					error('Could not prepare dataset for transfer:\n%s\n', ...
						output);
				end
				delete(cmd_file);
				
				[status, ~] = dos(['"' ppath ...
					'/tools/windows/pscp.exe" -batch -pw ' obj.Password ' '...
					obj.Username '@' obj.Hostname ':' cache_tmp ' ' ...
					pipeline_config.TempDir]);
				if status ~= 0, return, end
			end
			
			fprintf(1, [repmat('\b', 1, 13) repmat(' ', 1, 13) ...
				repmat('\b', 1, 13) '00%%']);
				
			cache_tmp = regexprep(cache_tmp, '.*/', '');
			
			untar([pipeline_config.TempDir '/' cache_tmp], ...
				pipeline_config.SampleCacheDir);
			delete([pipeline_config.TempDir '/' cache_tmp]);
		end
		
		function data = load_cached(obj, file)
			global pipeline_config;
			cache_path = [pipeline_config.SampleCacheDir '/' file];
			data = load(cache_path);
			delete(cache_path);
		end
		
		function data = load(obj, path)
			global pipeline_config;
			
			cache_path = [pipeline_config.SampleCacheDir '/' path];
			
			slashes = find(cache_path == '/');
			if isempty(slashes)
				cache_dir = '';
			else
				cache_dir = cache_path(1:slashes(end));
			end
			
			if ~exist(cache_dir)
				mkdir(cache_dir);
			end
			
			obj.prompt_password();
			
			if isunix || ismac
				[status, output] = unix(['scp -o BatchMode=true ' ...
					obj.Username '@' ...
					obj.Hostname ':' obj.Path '/' path ' ' cache_path]);
				if status ~= 0
					fprintf(1,'Could not load MAT file via SSH:\n');
					fprintf(1, '%s\n', output);
					return;
				end
			else
				[status, output] = dos(['"' ppath ...
					'/tools/windows/pscp.exe" -batch -pw ' obj.Password ' ' ...
					obj.Username '@' obj.Hostname ':' obj.Path '/' path ' ' ...
					cache_path]);
				if status ~= 0
					fprintf(1,'Could not load MAT file via SSH:\n');
					fprintf(1, '%s\n', output);
					return;
				end
			end
			
			data = load(cache_path);
		end
		
		function [] = prompt_password(obj)
			if ~ispc, return, end
			if isempty(obj.Password) || ~ischar(obj.Password)
				if 1
					obj.Password = passwordEntryDialog('WindowName', ...
						sprintf('SSH password for %s:', obj.Name), ...
						'CheckPasswordLength', false);
					if obj.Password == -1, obj.Password = ''; end
					if isempty(obj.Password)
						error(['You must specify an SSH password in ' ...
							'order to use an SSH repository.']);
					end
				else				
					obj.Password = input('Password: ', 's');
					ccount = length('Password: ') + length(obj.Password)+1;
					fprintf(1, [repmat('\b', 1, ccount) ...
						repmat(' ', 1, ccount) '\n']);
				end
			end
		end
	end
end

