function [] = export_json(json_file, varargin)

fid = fopen(json_file, 'W');
fprintf(fid, '{\n');

for k = 1:2:length(varargin)
	d = varargin{k+1};

	if numel(d) == 1
		fprintf(fid, '"%s": ', varargin{k});
		if isinteger(d) || islogical(d)
			fprintf(fid, '%d', d);
		elseif isreal(d)
			fprintf(fid, '%.3f', d);
		elseif ischar(d)
			fprintf(fid, '"%s"', d);
		elseif iscellstr(d)
			fprintf(fid, '["%s"]', d{1});
		else
			error 'Unsupported data type.';
		end
		if k < length(varargin)-1, fprintf(fid, ','); end
		fprintf(fid, '\n');
		continue;
	end
	
	fprintf(fid, '"%s": [', varargin{k});
			
	if length(d) ~= numel(d)
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
	else
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
	end
	fprintf(fid, ']');
	if k < length(varargin)-1, fprintf(fid, ','); end
	fprintf(fid, '\n');
end

fprintf(fid, '}\n');
fclose(fid);

