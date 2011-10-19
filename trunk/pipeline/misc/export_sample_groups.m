function [] = export_sample_groups(path, varargin)

[~, dir] = path_strip_dir(path);
[~, ~] = mkdir(dir);

fid = fopen(path, 'W');
fprintf(fid, '{\n');

for k = 1:2:length(varargin)
	if ~ischar(varargin{k}), error 'Invalid parameters.'; end
	
	fprintf(fid, '"%s": [\n', varargin{k});
	
	sample_ids = varargin{k+1};
	if isfield(sample_ids, 'Sample') && isfield(sample_ids.Sample, 'ID')
		sample_ids = sample_ids.Sample.ID;
	elseif isfield(sample_ids, 'Meta')
		sample_ids = sample_ids.Meta.Sample.ID;
	end
	
	if isempty(sample_ids)
		error('Empty sample group "%s".', varargin{k});
	end
	
	for s = 1:length(sample_ids)-1
		fprintf(fid, '\t"%s",\n', sample_ids{s});
	end
	fprintf(fid, '\t"%s"\n', sample_ids{end});
	
	if k < length(varargin)-1
		fprintf(fid, '],\n');
	else
		fprintf(fid, ']\n');
	end
end

fprintf(fid, '}\n');
fclose(fid);

