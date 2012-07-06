function order = aberration_grid(aberrations, image_file, color_method)

barh = 15;
gaph = 5;

if nargin < 3, color_method = 'red-gray-blue'; end

[F, S] = size(aberrations);

im = ones(F * (barh+gaph), S, 3);

for f = 1:F
	if isnumeric(color_method)
		val = aberrations(f, :);
		pixels = (1:barh) + (f-1) * (barh+gaph);
		for k = 1:max(val)
			im(pixels, val == k, 1) = color_method(k,1);
			im(pixels, val == k, 2) = color_method(k,2);
			im(pixels, val == k, 3) = color_method(k,3);
		end
		im(pixels, isnan(val), :) = 1;
		imwrite(im, image_file);
	end
end

