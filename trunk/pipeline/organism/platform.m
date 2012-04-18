function ret = platform(name, field)

global organism;

orgpath = [ppath '/organisms/' flatten_str(organism.Name) '/' ...
	flatten_str(organism.Version)];

if ~exist([orgpath '/platforms.mat'])
	error 'Platform descriptions were not found.';
end

platforms = load([orgpath '/platforms.mat']);
platforms = platforms.platforms;

platform = [];

for k = 1:length(platforms.Arrays)
	if regexpi(name, platforms.Arrays(k).Regex)
		platform = platforms.Arrays(k);
		break;
	end
end

if isempty(platform), error 'No platform by that name was found.'; end
	
if nargin == 1, ret = platform;	return; end

if strcmpi(field, 'probes')
	probes = load([ppath '/platforms/' platform.ProbeFile]);
	ret = probes.probes;
elseif strcmpi(field, 'cgh_probesets')
	probesets = load([ppath '/platforms/' platform.DefaultCghProbesets]);
	ret = probesets.probesets;
else
	error 'Unsupported platform field requested.';
end

