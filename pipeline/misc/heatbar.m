function [] = heatbar(data, color_method)

figure('PaperType', 'a4');
set(gca, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
data = data(:)';

if nargin == 1, color_method = 'redgreen'; end

if regexpi(color_method, 'primary')
	im = zeros(1, length(data), 3);
	im(1, data == 0 | data == 3, 1) = 1;
	im(1, data == 1 | data == 3, 2) = 1;
	im(1, data == 2, 3) = 1;
	image(im);
else
	if regexpi(color_method, 'redgreen')
		colormap(redgreencmap(256));
	elseif regexpi(color_method, 'gray')
		colormap(gray);
	elseif regexpi(color_method, 'hsv')
		colormap(hsv);
	end
	imagesc(data);
end

set(gca, 'Visible', 'off');


