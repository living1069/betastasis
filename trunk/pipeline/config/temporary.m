function path = temporary(base)

global pipeline_config;

% Note that tempname() returns paths like /tmp/fh9823fh289.
path = [pipeline_config.TempDir '/' base '_' datestr(now, 'mmmdd_HH:MM:SS') '/'];
[~,~] = mkdir(path);

