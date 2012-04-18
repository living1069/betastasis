function [] = pretty_bar(labels, values, varargin)

if isnumeric(labels)
	bar(labels, values, 'k');
elseif iscellstr(labels)
	bar(values, 0.5, 'k');
	set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels);
	xlim([0.5, length(labels)+0.5]);
	if any(cellfun(@length, labels) > 2)
		rotateticklabel(gca, 45);
	end
end


% Hide Y-axis
%set(gca, 'YTick', []);
%set(gca, 'YColor', 'w');



