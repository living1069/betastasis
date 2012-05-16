function [] = heatbar(data, color_method)

figure('PaperType', 'a4');
set(gca, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
data = data(:)';

if nargin == 1, color_method = 'green-black-red'; end

if iscell(color_method)
	for k = 1:max(data)
		im(1, data == k, :) = repmat(color_method{k}, sum(data == k), 1);
	end
	im(1, isnan(data), :) = repmat([.9 .9 .9], sum(isnan(data)), 1);
	image(im);
else
	if regexpi(color_method, 'green-black-red')
		colormap(redgreencmap(256));
	elseif regexpi(color_method, 'green-yellow-red')
		colormap(hsv2rgb([(1/3:-1/3/255:0)', ones(256, 2)]));
	elseif regexpi(color_method, 'gray')
		colormap(gray);
	elseif regexpi(color_method, 'hsv')
		colormap(hsv);
	end
	imagesc(data);
end

set(gca, 'Visible', 'off');


