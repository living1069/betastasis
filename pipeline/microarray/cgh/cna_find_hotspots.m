
% CNA_FIND_HOTSPOTS    Automatic identification of genomic CNA hotspots

% Author: Matti Annala <matti.annala@tut.fi>

function hotspots = cna_find_hotspots(segments, cgh_probesets, varargin)

global organism;

min_probes = 30;
min_significance = 5;
min_scale_step = 5;

for k = 1:2:length(varargin)
	if strcmpi(varargin{k}, 'MinProbes')
		min_probes = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'MinSignificance')
		min_significance = varargin{k+1};
		continue;
	end
	
	if strcmpi(varargin{k}, 'MinScaleStep')
		min_scale_step = varargin{k+1};
		continue;
	end
	
	error('Unrecognized option "%s".', varargin{k});
end

hotspots = struct;
hotspots.Chromosome = [];
hotspots.Position = [];
hotspots.Significance = [];

%sig = cna_significance_track(segments, cgh_probesets, ptemp);
%save ~/temp_sig.mat sig;
load ~/temp_sig.mat;

orig_sig = sig.LogP;

% We need to remove centromeres and other unprobed regions from consideration,
% because otherwise they will prevent whole chromosome losses from showing up
% as hotspots.
valid = ~isnan(orig_sig);
sig = orig_sig(valid);

cgh_probesets.Chromosome = cgh_probesets.Chromosome(valid);
cgh_probesets.Offset = cgh_probesets.Offset(valid);

max_amp_level = floor(max(sig))
max_del_level = floor(-min(sig))

hotspot_mask = zeros(size(sig));

progress = Progress;

% First check for significant hotspots of amplification.
for amp_level = max_amp_level:-5:min_significance
	sig_amp = medfilt1(double(sig > amp_level), min_probes);
	
	run_starts = [1; find(sig_amp(2:end) ~= sig_amp(1:end-1)) + 1; ...
		length(sig_amp) + 1];
	for r = 1:length(run_starts)-1
		run_start = run_starts(r);
		if sig_amp(run_start) == false, continue, end
		run_end = run_starts(r+1) - 1;
		
		run_len = run_end - run_start + 1;
		
		% Check if the hotspot is large enough.
		if run_len < min_probes, continue, end
			
		% Check if we have already found a more significant overlapping hotspot.
		% However, if the new region is min_scale_step times larger, then we
		% will include the new larger hotspot as well.
		if run_len < min_scale_step * max(hotspot_mask(run_start:run_end))
			continue;
		end
		
		% Found a new hotspot, store it.
		hotspots.Chromosome(end+1, 1) = cgh_probesets.Chromosome(run_start);
		hotspots.Position(end+1, :) = [cgh_probesets.Offset(run_start) ...
			cgh_probesets.Offset(run_end)];
		hotspots.Significance(end+1, 1) = amp_level;
		hotspot_mask(run_start:run_end) = run_len;
	end
	
	progress.update(0.5 * (max_amp_level - amp_level) / ...
		(max_amp_level - min_significance));
end

hotspot_mask = zeros(size(sig));

% Then check for significant hotspots of deletion.
for del_level = max_del_level:-5:min_significance
	sig_del = medfilt1(double(-sig > del_level), min_probes);
	
	run_starts = [1; find(sig_del(2:end) ~= sig_del(1:end-1)) + 1; ...
		length(sig_del) + 1];
	for r = 1:length(run_starts)-1
		run_start = run_starts(r);
		if sig_del(run_start) == false, continue, end
		run_end = run_starts(r+1) - 1;
		
		run_len = run_end - run_start + 1;
		
		% Check if the hotspot is large enough.
		if run_len < min_probes, continue, end
			
		% Check if we have already found a more significant overlapping hotspot.
		% However, if the new region is min_scale_step times larger, then we
		% will include the new larger hotspot as well.
		if run_len < min_scale_step * max(hotspot_mask(run_start:run_end))
			continue;
		end

		% Found a new hotspot, store it.
		hotspots.Chromosome(end+1, 1) = cgh_probesets.Chromosome(run_start);
		hotspots.Position(end+1, :) = [cgh_probesets.Offset(run_start) ...
			cgh_probesets.Offset(run_end)];
		hotspots.Significance(end+1, 1) = -del_level;
		hotspot_mask(run_start:run_end) = run_len;
	end
	
	progress.update(0.5 + 0.5 * (max_del_level - del_level) / ...
		(max_del_level - min_significance));
end

cnv.Chromosome = hotspots.Chromosome;
cnv.CNVType = (hotspots.Significance > 0) + 1;
cnv.Start = hotspots.Position(:, 1);
cnv.End = hotspots.Position(:, 2);

hotspots

hs_cytobands = cytobandread('~/organisms/homo_sapiens/cytoBand.txt');
figure; chromosomeplot(hs_cytobands, 'CNV', cnv);
saveas(gcf, '~/hotspots.pdf');

