
% Author: Matti Annala <matti.annala@tut.fi>

function [] = tophat_analysis(reads, varargin)

global organism;
exons = organism.Exons;

seq_files = seq_resource_files(reads);

junctions = struct;
junctions.Junctions = cell(1, length(seq_files));

%for k = 1:length(seq_files)
for k = 1
	tophat_out_dir = ptemp;
	
	[flags, index_suffix] = tophat_flags_for_reads(seq_files{k});
	
	% FIXME: Enable --allow-indels?
	status = unix(sprintf([ ...
		'PATH=' ppath '/tools/bowtie:' ppath '/tools/samtools:$PATH ' ...
		'%s/tools/tophat/tophat %s -o %s -a8 %s%s %s '], ...
		ppath, flags, tophat_out_dir, ...
		bowtie_index('genome'), index_suffix, seq_files{k}));
	if status ~= 0, error 'Tophat junction discovery failed.'; end
	
	junctions.Junctions{k} = struct;
	
	fid = fopen([tophat_out_dir '/junctions.bed']);
	data = textscan(fid, '%s %*d %*d %*s %d %*s %*s %*s %*s %d %s %s');
	fclose(fid);
	
	chromosomes = chromosome_sym2num(data{1});
	scores = data{2};
	block_counts = data{3};
	block_sizes = data{4};
	block_starts = data{5};
	
	
end










function [flags, index_suffix] = tophat_flags_for_reads(reads)

[color, quality] = seq_read_type(reads);

flags = '';
index_suffix = '';
if color == 1, index_suffix = '_colorspace'; end

if color == 1 && quality == 1
	flags = '-C --solexa-quals --integer-quals';
elseif color == 1 && quality == 0
	flags = '-C';
elseif color == 0 && quality == 1
	flags = '';
elseif color == 0 && quality == 0
	flags = '';
end

