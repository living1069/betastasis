function [] = aberration_grid(image_file, aberrations, varargin)

order = 1:size(aberrations{1,1}, 2);
box_size = [14 14];

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'order')
		order = varargin{k+1}; continue;
	end
	if rx(varargin{k}, 'size')
		box_size = varargin{k+1}; continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

colors = color(aberrations{1,1}, aberrations{1,2});
for k = 2:size(aberrations, 1)
	colors = cat(1, colors, color(aberrations{k, 1}, aberrations{k, 2}));
end

colors = colors(:, order);

[F, S] = size(colors); F = F / 3;

boxw = box_size(1);
boxh = box_size(2);
gapw = 1;
gaph = 1;

%if rx(colors, 'no.*h.*gap')
%	boxw = 1;
%	gapw = 0;
%end
	
% Initialize the bitmap with a white background...
im = ones(F * (boxh+gaph) + gaph, S * (boxw+gapw) + 2 - gapw, 3);

% And a black border.
im(1, :, :) = 0;
im(end, :, :) = 0;
im(:, 1, :) = 0;
im(:, end, :) = 0;

for f = 1:F
	ypixels = (1:boxh) + (f-1) * (boxh+gaph) + gaph;
	
	for s = 1:S
		xpixels = (1:boxw) + (s-1) * (boxw+gapw) + 1;
		im(ypixels, xpixels, 1) = colors(f*3-2, s);
		im(ypixels, xpixels, 2) = colors(f*3-1, s);
		im(ypixels, xpixels, 3) = colors(f*3, s);
	end
	
	imwrite(im, image_file);
end














function colors = color(data, palette)

colors = zeros(size(data, 1) * 3, size(data, 2));

missing = [0 0 1];
neutral = [0 0 .9];

if rx(palette, 'blue-red')
	if any(abs(data(:)) > 2), error('Data must be between -2 and 2.'); end
	palette = hsv2rgb([missing; 2/3 1 1; 2/3 .5 1; neutral; 0 .5 1; 0 1 1]);
	data = data + 4;
	
elseif rx(palette, 'green')
	if any(data(:) > 2 | data(:) < 0)
		error('Data must be between 0 and 2.');
	end
	palette = hsv2rgb([missing; neutral; 1/3 .5 1; 1/3 1 1]);
	data = data + 2;
	
elseif rx(palette, 'yellow')
	if any(data(:) > 2 | data(:) < 0)
		error('Data must be between 0 and 2.');
	end
	palette = hsv2rgb([missing; neutral; 1/6 .5 1; 1/6 1 1]);
	data = data + 2;
	
elseif rx(palette, 'black')	
	if any(data(:) > 2 | data(:) < 0)
		error('Data must be between 0 and 2.');
	end
	palette = hsv2rgb([missing; neutral; 0 0 .4; 0 0 0]);
	data = data + 2;

elseif rx(palette, 'rgby')
	if any(data(:) > 4 | data(:) < 0)
		error('Data must be between 0 and 4.');
	end
	palette = hsv2rgb([missing; neutral; ...
		[0 128 230; 80 128 230; 150 128 230; 40 140 250]/255]);
	data = data + 2;

end

data(isnan(data)) = 1;    % Index 1 is reserved for missing values.

for k = 1:size(data, 1)
	colors(k*3-2, :) = palette(data(k, :), 1);
	colors(k*3-1, :) = palette(data(k, :), 2);
	colors(k*3, :)   = palette(data(k, :), 3);
end

