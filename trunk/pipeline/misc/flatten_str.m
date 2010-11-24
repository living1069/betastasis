function str = flatten_str(str)
str = strrep(str, ' ', '_');
str = strrep(str, '.', '_');
str = strrep(str, '-', '_');
str = lower(str);

