function [] = gene_survival(expr, gene, threshold, file_prefix, varargin)

criterion = 'survival';

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'Criterion')
		criterion = varargin{k+1};
		continue;
	end

	error('Unrecognized option "%s".', varargin{k});
end

idx = gene_idx(gene);

log_expr = log2(expr.Mean(idx, :));

if regexpi(criterion, 'survival')
	if ~isfield(expr.Meta.Patient, 'SurvivalTime')
		error 'Cannot find survival time information.';
	end
	survival = expr.Meta.Patient.SurvivalTime;
	censored = expr.Meta.Patient.Censored;
	label = 'Percent survival';
elseif regexpi(criterion, 'recurrence')
	if ~isfield(expr.Meta.Patient, 'RecurrenceFreeTime')
		error 'Cannot find recurrence time information.';
	end
	survival = expr.Meta.Patient.RecurrenceFreeTime;
	censored = ~expr.Meta.Patient.RecurrenceEvent;
	label = 'Percent recurrence free';
else
	error 'Unsupported criterion.';
end

valid = ~isnan(survival) & ~isnan(censored);
log_expr = log_expr(valid);
survival = survival(valid);
censored = censored(valid);

time_label = 'Days';
if max(survival) > 2000
	survival = survival / 30;
	time_label = 'Months';
end

low = (log_expr < threshold);
high = (log_expr >= threshold);

figure('PaperPosition', [.1 .1 .8 .85], 'PaperUnits', 'normalized');
subplot('position', [.2 .05 .6 .15]); hold all;
hist(log_expr, 100);
h = findobj(gca, 'Type', 'patch'); set(h, 'FaceColor', [.8 .8 .8]);
xlabel([gene ' expression (log2)'], 'Interpreter', 'none');
ylabel('Number of samples');

plot([threshold threshold], [0 max(ylim)], 'r', 'LineWidth', 2);

subplot('position', [.1 .25 .8 0.7]);
p = logrank(survival(low), censored(low), survival(high), censored(high));

if p >= 0.0001
	p_str = sprintf('%.5f', p);
else
	p_str = sprintf('%.3e', p);
end

title(sprintf('Kaplan-Meier survival plot for %s (p-value %s, logrank test)',...
	gene, p_str), 'Interpreter', 'none');
xlabel([time_label ' since initial diagnosis']); ylabel(label);
legend('Low expression', 'High expression', 'Censored');


saveas(gcf, file_prefix);

