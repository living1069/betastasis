function [] = survival_plot(survtimes, censored)

if nargin == 1
	if isfield(survtimes, 'Meta') && isfield(survtimes.Meta, 'Patient')
		patient = survtimes.Meta.Patient;
	elseif isfield(survtimes, 'Patient')
		patient = survtimes.Patient;
	else
		error 'No survival time metadata available.';
	end
	
	survtimes = patient.SurvivalTime;
	censored = patient.Censored;
end

[f, x, flo, fup] = ecdf(survtimes, 'censoring', censored, ...
	'function', 'survivor');

hold all;
stairs(x, f);
%stairs(x, flo, 'r:');
%stairs(x, fup, 'r:');
xlabel('Days since initial diagnosis');
ylabel('Probability of survival');

