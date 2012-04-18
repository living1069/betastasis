function [] = correlate_with_clinical(dataset, pos_samples, neg_samples, labels)

global organism;
genes = organism.Genes;
transcripts = organism.Transcripts;
exons = organism.Exons;

meta = dataset;
if isfield(meta, 'Meta'), meta = meta.Meta; end

if iscellstr(pos_samples)
	pos_samples = ismember(meta.Sample.ID, pos_samples);
end

if nargin == 2
	neg_samples = ~pos_samples;
elseif iscellstr(neg_samples)
	neg_samples = ismember(meta.Sample.ID, neg_samples);
end

if nargin < 4
	labels = { 'Positive', 'Negative' };
end

cecum = strcmpi('cecum', meta.Misc.anatomic_organ_subdivision);
ascending_colon = strcmpi('ascending colon', ...
	meta.Misc.anatomic_organ_subdivision);
hepatic_flexure = strcmpi('hepatic flexure', ...
	meta.Misc.anatomic_organ_subdivision);
transverse_colon = strcmpi('transverse colon', ...
	meta.Misc.anatomic_organ_subdivision);
descending_colon = strcmpi('descending colon', ...
	meta.Misc.anatomic_organ_subdivision);
sigmoid_colon = strcmpi('sigmoid colon', meta.Misc.anatomic_organ_subdivision);
rectum = strcmpi('rectum', meta.Misc.anatomic_organ_subdivision);

stage_i = rx(meta.Misc.tumor_stage, 'Stage I[^I]*$');
stage_ii = rx(meta.Misc.tumor_stage, 'Stage II[^I]*$');
stage_iii = rx(meta.Misc.tumor_stage, 'Stage III[^I]*$');
stage_iv = strcmpi('M1', meta.Misc.distant_metastasis_pathologic_spread);
stage_known = stage_i | stage_ii | stage_iii | stage_iv;

msi = strcmpi(meta.Misc.microsatellite_instability, 'YES');
mss = strcmpi(meta.Misc.microsatellite_instability, 'NO');

lymph_invasion = strcmpi(meta.Misc.lymphatic_invasion_present, 'YES');
no_lymph_invasion = strcmpi(meta.Misc.lymphatic_invasion_present, 'NO');

vascular_invasion = strcmpi('YES', meta.Misc.vascular_invasion_present);
no_vascular_invasion = strcmpi('NO', meta.Misc.vascular_invasion_present);

polyp_history = strcmpi(meta.Misc.history_of_colon_polyps, 'YES');
no_polyp_history = strcmpi(meta.Misc.history_of_colon_polyps, 'NO');

herited = rx(meta.Misc.number_of_first_degree_relatives_with_cancer_diagnosis, '^[^0-]$');
not_herited = rx(meta.Misc.number_of_first_degree_relatives_with_cancer_diagnosis, '^0$');

male = strcmpi(meta.Patient.Gender, 'male');
female = strcmpi(meta.Patient.Gender, 'female');

plot_freq_in_subgroups({stage_i, stage_ii, stage_iii, stage_iv, ...
	msi, mss, lymph_invasion, no_lymph_invasion, ...
	vascular_invasion, no_vascular_invasion, ...
	polyp_history, no_polyp_history, ...
	cecum, ascending_colon, hepatic_flexure, transverse_colon, ...
	descending_colon, sigmoid_colon, rectum, ...
	herited, not_herited, male, female }, ...
	{ 'Stage I', 'Stage II', 'Stage III', 'Stage IV', ...
	'MSI', 'MSS', ...
	'Lymph invasion', 'No lymph invasion', ...
	'Vascular invasion', 'No vascular invasion', ...
	'Polyp history', 'No polyp history', ...
	'Cecum', 'Ascending colon', 'Hepatic flexure', 'Transverse colon', ...
	'Descending colon', 'Sigmoid colon', 'Rectum', ...
	'Herited', 'Not herited', 'Male', 'Female' }, ...
	'Fusion positive', pos_samples, 'Fusion negative', neg_samples);


