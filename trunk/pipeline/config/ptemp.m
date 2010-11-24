function p = ptemp()

global pipeline_config;

% Note that tempname() returns paths like /tmp/fh9823fh289.
if isunix || ismac
	p = [pipeline_config.TempDir strrep(tempname, '/tmp', '')];
else
	p = [pipeline_config.TempDir regexprep(tempname, '^(.+)(\\.+?)$', '$2')];
end

