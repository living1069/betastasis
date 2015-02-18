function [] = heatmap(im, filename, varargin)

cmap_name = 'blue.*red';
domain = [];

for k = 1:2:length(varargin)
	if rx(varargin{k}, 'domain')
		domain = varargin{k+1}; continue;
	end
	if rx(varargin{k}, 'color')
		cmap_name = varargin{k+1}; continue;
	end
	error('Unrecognized option "%s".', varargin{k});
end

if isempty(domain), domain = [nanmin(im(:)), nanmax(im(:))]; end
if length(domain) == 2
	domain = [domain(1), (domain(1)+domain(2))/2, domain(2)];
end

colormap = zeros(256, 3);
if rx(cmap_name, 'blue.*red')
	colormap(1:128, 1) = 2/3; colormap(129:256, 1) = 0;
	colormap(1:128, 2) = linspace(1, 0, 128);
	colormap(129:256, 2) = linspace(0, 1, 128);
	colormap(1:128, 3) = linspace(1, .9, 128);
	colormap(129:256, 3) = linspace(.9, 1, 128);
	colormap = hsv2rgb(colormap);
elseif rx(cmap_name, 'gray.*red')
	colormap(:, 1) = 0;
	colormap(:, 2) = linspace(0, 1, 256);
	colormap(:, 3) = 0.9;
	colormap = hsv2rgb(colormap);
else
	error('Unrecognized colormap requested.');
end

im(im < domain(1)) = domain(1);
im(im > domain(end)) = domain(end);
im = interp1(domain, [1 128 256], im, 'linear', 'extrap');
imwrite(im, colormap, filename);


