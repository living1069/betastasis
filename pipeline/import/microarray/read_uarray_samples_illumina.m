function [raw, probesets] = read_uarray_samples_illumina(sample_file)

global organism;

[data, headers] = readtable(sample_file, 'Numeric', {'AVG_Signal'});

probe_seq = data{rx(headers, 'PROBE_SEQUENCE')};

S = sum(rx(headers, 'AVG_Signal'));
fprintf(1, 'Found %d samples.\n', S);

probes.Sequence = probe_seq;
probesets = create_gene_probesets(probes, 'EveryTranscript', false, ...
	'MinProbesetSize', 0);
	
raw.mean = cat(2, data{rx(headers, 'AVG_Signal')});
raw.meta.type = 'Microarray probe intensities';
raw.meta.platform = 'Illumina';
raw.meta.sample_id = regexprep(headers(rx(headers, 'AVG_Signal')), ...
	'(.+?).AVG_Signal', '$1');

