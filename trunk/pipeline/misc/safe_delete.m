function [] = safe_delete(file)

global pipeline_config;

if isempty(strfind(file, pipeline_config.TempDir))
	fprintf(1, ['WARNING: Omitted deleting file %s that is stored outside ' ...
	            'the temporary data directory.\n']);
	return;
end

old_state = recycle('off');
delete(file);
recycle(old_state);

