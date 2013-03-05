function order = aberration_grid(aberrations, image_file, colors)

if nargin < 3, colors = 'blue-red'; end

[F, S] = size(aberrations);

boxw = round((400 - S) / S);
boxh = 14;
gapw = 1;
gaph = 1;

if rx(colors, 'no.*h.*gap')
	boxw = 1;
	gapw = 0;
end
	
if rx(colors, 'blue.*red')
	% 1..5 = blue to red discretized scale
	% 6..7 = light green, dark green (mutations)
	colors = hsv2rgb([2 3 3; 2 1.5 3; 2 0 2.7; 0 1.5 3; 0 3 3; ...
		1 1.5 3; 1 3 3]/3);
end

% Initialize the bitmap with a white background...
im = ones(F * (boxh+gaph) + gaph, S * (boxw+gapw) + 2 - gapw, 3);

% And a black border.
im(1, :, :) = 0;
im(end, :, :) = 0;
im(:, 1, :) = 0;
im(:, end, :) = 0;

for f = 1:F
	val = aberrations(f, :);
		
	ypixels = (1:boxh) + (f-1) * (boxh+gaph) + gaph;
	
	for s = 1:S
		xpixels = (1:boxw) + (s-1) * (boxw+gapw) + 1;
		if isnan(val(s))
			im(ypixels, xpixels, :) = 1;
		else
			im(ypixels, xpixels, 1) = colors(val(s), 1);
			im(ypixels, xpixels, 2) = colors(val(s), 2);
			im(ypixels, xpixels, 3) = colors(val(s), 3);
		end
	end
	
	imwrite(im, image_file);
end

