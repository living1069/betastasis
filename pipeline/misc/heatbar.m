function [] = heatbar(data, image_file, color_method)

barh = 15;

data = data(:)';

if nargin < 3, color_method = 'green-black-red'; end
	
im = ones(barh, length(data), 3);

if isnumeric(color_method)
	for k = 1:max(data)
		im(:, data == k, 1) = color_method(k,1);
		im(:, data == k, 2) = color_method(k,2);
		im(:, data == k, 3) = color_method(k,3);
	end
	im(1, isnan(data), :) = repmat([.9 .9 .9], sum(isnan(data)), 1);
else
	if regexpi(color_method, 'green-black-red')
		cmap = redgreencmap(256);
	elseif regexpi(color_method, 'green-yellow-red')
		cmap = hsv2rgb([(1/3:-1/3/255:0)', ones(256, 2)]);
	elseif regexpi(color_method, 'gray')
		cmap = gray;
	elseif regexpi(color_method, 'hsv')
		cmap = hsv;
	else
		error 'Unknown colormap specified.';
	end
	
	scaled = data - min(data);
	scaled = scaled / max(scaled) * 255;
	im(:, :, 1) = repmat(cmap(round(scaled+1), 1)', barh, 1);
	im(:, :, 2) = repmat(cmap(round(scaled+1), 2)', barh, 1);
	im(:, :, 3) = repmat(cmap(round(scaled+1), 3)', barh, 1);
end

imwrite(im, image_file);

