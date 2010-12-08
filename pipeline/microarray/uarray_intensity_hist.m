function [] = uarray_intensity_hist(samples, image_prefix, varargin)

if isstruct(samples), samples = samples.Mean; end
	
log_transform = true;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'LogTransform')
		log_transform = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

S = size(samples, 2);

for s = 1:S
	if any(samples(:, s) <= 0)
		fprintf(1, 'WARNING: Sample #%d contains %d non-positive values.\n', ...
			s, sum(samples(:, s) <= 0));
	end
end

if log_transform
	samples = log2(samples);
end

ninf = (samples == -Inf);
samples(ninf) = NaN;
limits = [floor(min(min(samples))) ceil(max(max(samples)))];
samples(ninf) = -Inf;

bins = limits(1):(limits(2) - limits(1))/500:limits(2);

for s = 1:S
	figure; hist(samples(:, s), bins);
	xlim(limits); ylabel('Number of probes');
	
	if log_transform
		xlabel('Probe log-intensity');
	else
		xlabel('Probe intensity');
	end
	
	if S == 1
		saveas(gcf, [image_prefix '.pdf']);
	else
		saveas(gcf, [image_prefix '_' num2str(s) '.pdf']);
	end
end


