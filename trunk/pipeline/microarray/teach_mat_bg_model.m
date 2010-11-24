function background = teach_mat_bg_model(samples, probes)

fprintf(1, 'Calculating MAT background model parameters...\n');

orig_probes = probes;

% Pick the lowest 20% of probes by intensity for teaching the model. The
% underlying assumption is that these low intensity probes are only measuring
% non-specific hybridization (i.e. the background).
bg_probes = find(median(samples.Mean, 2) < ...
	quantile(median(samples.Mean, 2), 0.2));
probes = struct('Sequence', { probes.Sequence(bg_probes) });

H = mat_model_probe_matrix(probes);

% The least squares solution is calculated blockwise because the observation
% matrix H consists of vertically repeated blocks.
% Tikhonov regularization is used because H is rank deficient.
x = samples.Mean(bg_probes, :);
theta = inv(H' * H + 0.1 * eye(size(H, 2))) * H' * mean(log2(x), 2);

% This is numerically unstable because H is very rank deficient.
%theta = H \ mean(log2(x), 2);

% Render some plots of the model parameters.
beta_a = theta(1:25);
beta_c = theta(26:50);
beta_g = theta(51:75);
beta_t = theta(76:100);

figure; hold all;
plot(beta_a);
plot(beta_c);
plot(beta_g);
plot(beta_t);
legend('Beta-A', 'Beta-C', 'Beta-G', 'Beta-T');
saveas(gcf, '~/mat_beta.png');

fprintf(1, 'Gamma-A coefficient: %f\n', theta(101));
fprintf(1, 'Gamma-C coefficient: %f\n', theta(102));
fprintf(1, 'Gamma-G coefficient: %f\n', theta(103));
fprintf(1, 'Gamma-T coefficient: %f\n', theta(104));

background = mat_model_probe_matrix(orig_probes) * theta;

