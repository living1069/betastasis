function [] = export_json(json_files, varargin)

if ischar(json_files) && size(json_files, 1) == 1
	json_files = { json_files };
end

N = length(json_files);

for f = 1:N
	fid = fopen(json_files{f}, 'W');
	
	fprintf(fid, '{\n');
	
	for k = 1:2:length(varargin)
		d = varargin{k+1};

		if size(d, 2) == N && size(d, 1) == 1
			fprintf(fid, '"%s": ', varargin{k});
			if isinteger(d) || islogical(d)
				fprintf(fid, '%d', d(f));
			elseif isreal(d)
				fprintf(fid, '%.3f', d(f));
			elseif iscellstr(d) 
				fprintf(fid, '"%s"', d{f});
			end
			if k < length(varargin)-1, fprintf(fid, ','); end
			fprintf(fid, '\n');
			continue;
		end
		
		fprintf(fid, '"%s": [', varargin{k});
				
		if N == 1 && length(d) ~= numel(d)
			% The field is to be stored as a JSON array of arrays (matrix).
			for r = 1:size(d, 1)
				fprintf(fid, '[');
				if isinteger(d) || islogical(d)
					fprintf(fid, '%d, ', d(r, 1:end-1));
					fprintf(fid, '%d', d(r, end));
				elseif isreal(d)
					fprintf(fid, '%.3f, ', d(r, 1:end-1));
					fprintf(fid, '%.3f', d(r, end));
				elseif iscellstr(d)
					for c = 1:size(d, 2)-1
						fprintf(fid, '"%s", ', d{r, c});
					end
					fprintf(fid, '"%s"', d{r, end});
				end
				fprintf(fid, ']');
				if r < size(d, 1), fprintf(fid, ', '); end
			end
		elseif N == 1
			% The field is to be stored as a JSON array.
			if isinteger(d) || islogical(d)
				if length(d) > 1
					fprintf(fid, '%d, ', d(1:end-1));
				end
				fprintf(fid, '%d', d(end));
			elseif isreal(d)
				if length(d) > 1
					fprintf(fid, '%.3f, ', d(1:end-1));
				end
				fprintf(fid, '%.3f', d(end));
			elseif iscellstr(d)
				for c = 1:length(d)-1
					fprintf(fid, '"%s", ', d{c});
				end
				fprintf(fid, '"%s"', d{end});
			end
		elseif size(d, 2) == N
			% One slice of the field is to be stored as a JSON array.
			if isinteger(d) || islogical(d)
				fprintf(fid, '%d, ', d(1:end-1, f));
				fprintf(fid, '%d', d(end, f));
			elseif isreal(d)
				fprintf(fid, '%.3f, ', d(1:end-1, f));
				fprintf(fid, '%.3f', d(end, f));
			elseif iscellstr(d) 
				for c = 1:N-1
					fprintf(fid, '"%s", ', d{c, f});
				end
				fprintf(fid, '"%s"', d{end, f});
			end
		end
		fprintf(fid, ']');
		if k < length(varargin)-1, fprintf(fid, ','); end
		fprintf(fid, '\n');
	end

	fprintf(fid, '}\n');
	fclose(fid);
end


