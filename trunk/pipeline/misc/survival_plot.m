function [] = survival_plot(survival, censored)

if nargin == 1
	if isfield(survival, 'meta') && isfield(survival.meta, 'survival_time')
		censored = survival.meta.survival_time_censored;
		survival = survival.meta.survival_time;
	elseif isfield(survival, 'survival_time')
		censored = survival.survival_time_censored;
		survival = survival.survival_time;
	else
		error 'No survival time metadata available.';
	end
end

[f, x, flo, fup] = ecdf(survival, 'censoring', censored, ...
	'function', 'survivor');

hold all;
stairs(x, f);
%stairs(x, flo, 'r:');
%stairs(x, fup, 'r:');
xlabel('Days since initial diagnosis');
ylabel('Probability of survival');

