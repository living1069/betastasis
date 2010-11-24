
% RENDER_UARRAY_INTENSITIES     Renders a spatial map of probe intensities
%
%    RENDER_UARRAY_INTENSITIES(SAMPLES, PROBES, IMAGE_PREFIX) renders PNG
%    images of the raw microarray probe intensity data in SAMPLES. The probe
%    configuration of the microarray must also be provided as the second
%    argument PROBES. The PNG images will be written to disk using IMAGE_PREFIX,
%    so that if IMAGE_PREFIX = 'figure', the images will have filenames
%    'figure_1.png', 'figure_2.png', etc.
%
%    Intensities will be rendered using the natural scale. Probes with unknown
%    intensity or NaN intensity will be rendered as black pixels.
%    
%    RENDER_UARRAY_INTENSITIES(..., 'RedMissing', true) tells the function to
%    show spots with unknown or NaN intensity as red pixels.

function [] = render_uarray_intensities(samples, probes, image_prefix, varargin)

red_missing = false;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'RedMissing')
		red_missing = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

if isstruct(samples)
	samples = samples.Mean;
end

window_size = 7;
window_elems = window_size^2;
radius = (window_size - 1) / 2;

probes.XPos = probes.XPos + (1 - min(probes.XPos));
probes.YPos = probes.YPos + (1 - min(probes.YPos));

h = max(probes.YPos);
w = max(probes.XPos);

progress = Progress;

for s = 1:size(samples, 2)
	f = nan(h + 2*radius, w + 2*radius);
	for p = 1:size(samples, 1)
		f(probes.YPos(p) + radius, probes.XPos(p) + radius) = samples(p, s);
	end

	ignore = find(isnan(f));
	f(ignore) = (mod(ignore, 2) - 0.5) * Inf;

	f = medfilt2(f, [window_size window_size]);
	f(ignore) = NaN;

	f = f(radius+1:radius+h, radius+1:radius+w);
	f(f == Inf) = NaN;
	f(f == -Inf) = NaN;
		
	image = f - min(min(f));
	image = image / max(max(image));
	image = repmat(image, [1 1 3]);
	napix = find(isnan(image(:, :, 1)));
	plane_size = size(image, 1) * size(image, 2);
	image(napix) = red_missing;      % Black or red depending on user choice.
	image(napix + plane_size) = 0;
	image(napix + 2*plane_size) = 0;
	
	imwrite(image, [image_prefix '_' num2str(s) '.png']);
	progress.update(s / size(samples, 2));
end

