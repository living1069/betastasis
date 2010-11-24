function H = mat_model_probe_matrix(probes)

N = length(probes.Sequence);

% It is assumed that all probes are of equal length. This restriction could
% possibly be relaxed by improving the model.
probe_len = length(probes.Sequence{1});
for k = 2:N
	if length(probes.Sequence{k}) ~= probe_len
		error 'Probes must be of equal length.';
	end
end

% Indicator function of whether a nucleotide at some position in a sequence
% is of a certain type.
I_jk = zeros(N, 4 * probe_len);
for p = 1:N
	I_jk(p, 1:25) = probes.Sequence{p} == 'a';
	I_jk(p, 26:50) = probes.Sequence{p} == 'c';
	I_jk(p, 51:75) = probes.Sequence{p} == 'g';
	I_jk(p, 76:100) = probes.Sequence{p} == 't';
end

% Squared nucleotide counts.
nsq = zeros(N, 4);
for p = 1:N
	nsq(p, 1) = sum(probes.Sequence{p} == 'a')^2;
	nsq(p, 2) = sum(probes.Sequence{p} == 'c')^2;
	nsq(p, 3) = sum(probes.Sequence{p} == 'g')^2;
	nsq(p, 4) = sum(probes.Sequence{p} == 't')^2;
end

% Combined observation matrix.
H = [ I_jk nsq ones(N, 1) ];

