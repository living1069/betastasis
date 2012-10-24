function order = aberration_grid(aberrations, image_file, color_method)

[F, S] = size(aberrations);

boxw = round((400 - S) / S);
boxh = 14;
gap = 1;

% Initialize the bitmap with a light gray background...
im = ones(F * (boxh+gap) + gap, S * (boxw+gap) + gap, 3) * 0.95;

% A black border...
for f = 1:F-1, im(gap + f * (boxh+gap), :, :) = 1; end
for s = 1:S-1, im(:, gap + s * (boxw+gap), :) = 1; end

% And a grid of white lines.
im(1, :, :) = 0;
im(end, :, :) = 0;
im(:, 1, :) = 0;
im(:, end, :) = 0;

for f = 1:F
	val = aberrations(f, :);
		
	ypixels = (1:boxh) + (f-1) * (boxh+gap) + gap;
	
	for s = 1:S
		xpixels = (1:boxw) + (s-1) * (boxw+gap) + gap;
		if isnan(val(s))
			im(ypixels, xpixels, :) = 1;
		else
			im(ypixels, xpixels, 1) = color_method(val(s), 1);
			im(ypixels, xpixels, 2) = color_method(val(s), 2);
			im(ypixels, xpixels, 3) = color_method(val(s), 3);
		end
	end
	
	imwrite(im, image_file);
end

