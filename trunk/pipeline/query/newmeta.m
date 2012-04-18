
function b = newmeta(a)

if isfield(a, 'meta')
	error 'Looks like the dataset is already in new format.';
end

b = struct;
if isfield(a, 'Mean'), b.mean = a.Mean; end
if isfield(a, 'Fusions'), b.fusions = a.Fusions; end

if isfield(a, 'Meta')
	if isfield(a.Meta, 'Type')
		b.meta.type = a.Meta.Type;
	end
	
	if isfield(a.Meta, 'Sample')
		if isfield(a.Meta.Sample, 'ID')
			b.meta.sample_id = a.Meta.Sample.ID';
		end
		if isfield(a.Meta.Sample, 'Filename')
			b.meta.uarray_filename = a.Meta.Sample.Filename';
		end
		if isfield(a.Meta.Sample, 'Channel')
			b.meta.uarray_channel = a.Meta.Sample.Channel';
		end
		if isfield(a.Meta.Sample, 'Type')
			b.meta.sample_type = a.Meta.Sample.Type';
		end
		if isfield(a.Meta.Sample, 'Organ')
			b.meta.organ = a.Meta.Sample.Organ';
		end
	end
	
	if isfield(a.Meta, 'Patient')
		if isfield(a.Meta.Patient, 'ID')
			b.meta.patient_id = a.Meta.Patient.ID';
		end
		if isfield(a.Meta.Patient, 'Gender')
			b.meta.gender = a.Meta.Patient.Gender';
		end
		if isfield(a.Meta.Patient, 'SurvivalTime')
			b.meta.survival_time = a.Meta.Patient.SurvivalTime';
		end
		if isfield(a.Meta.Patient, 'Censored')
			b.meta.survival_time_censored = a.Meta.Patient.Censored';
		end
		if isfield(a.Meta.Patient, 'Status')
			b.meta.patient_status = a.Meta.Patient.Status';
		end
	end
end
