
% This function reads raw probe intensities from a .CEL file. The mean
% and standard deviation is returned for every known probe in the microarray.
% These values are originally estimated by the microarray platform itself.
% 
% Inputs:
%     cel_file - File path to the .CEL file that is to be read.
%     probes - The set of probes whose intensities are to be read.

% Author: Matti Annala <matti.annala@tut.fi>

function sample = read_uarray_sample_cel(cel_file, probes)

N = length(probes.Sequence);

sample = struct;
sample.mean = zeros(N, 1, 'single');
	
cel = affyread(cel_file);

idx = probes.YPos * cel.Cols + probes.XPos + 1;
sample.mean = cel.Probes(idx, 3);
%sample.Stdev = cel.Probes(idx, 4);

