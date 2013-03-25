function files = find_files(pattern, root)

if nargin < 2, root = '.'; end

[status, out] = unix(sprintf('find %s -not -type d', root));
if status ~= 0, error 'UNIX find command failed.'; end

data = textscan(out, '%s', 'Delimiter', '\n');
files = data{1};

files = files(rx(files, pattern));
files = sort_nat(files);


