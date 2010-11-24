function sample = gmedian_normalize(sample, ref)

smed = median(sample);
sample = sample - smed;
sample = sample / std(sample) * std(ref);
sample = sample + median(ref);

