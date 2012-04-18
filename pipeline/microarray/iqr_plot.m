function outliers = iqr_plot(raw, image_file, varargin)

S = size(raw.Mean, 2);

order_by = 'Median';
marked = false(1, S);

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Order')
		order_by = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'Mark')
		if islogical(varargin{k+1})
			marked = varargin{k+1};
		elseif isnumeric(varargin{k+1})
			marked(varargin{k+1}) = true;
		else
			error 'Marked samples specified in an invalid manner.';
		end
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

qtiles = quantile(raw.Mean, [.05 .25 .5 .75 .95]);
sample_names = raw.Meta.Sample.ID;

if strcmpi(order_by, 'Median')
	[~, order] = sort(qtiles(3, :));
elseif strcmpi(order_by, 'IQR')
	[~, order] = sort(qtiles(4, :) - qtiles(2, :));
elseif isempty(order_by)
	order = 1:S;
else
	error 'Invalid ordering requested.';
end

marked = find(marked(order));

sample_names = sample_names(order);
qtiles = qtiles(:, order);

figure; hold all;
for k = 1:5, plot(qtiles(k, :), 'k', 'LineWidth', 2); end

scatter(marked, qtiles(3, marked), 32, 'r', 'filled');

if length(raw.Meta.Sample.ID) <= 10
	set(gca, 'xtick', 1:length(raw.Meta.Sample.ID), ...
		'xticklabel', raw.Meta.Sample.ID);
	rotateticklabel(gca, 45);
	xlabel('Sample ID');
else
	set(gca, 'xtick', []);
end

ylabel('Quantile traces');

ylim([0, median(qtiles(5, :)) * 2]);

saveas(gcf, image_file);

