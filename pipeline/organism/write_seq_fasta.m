
% This function writes out all transcript transcripts in a transcriptome into 
% a FASTA format file. Each sequence will be preceded by a header line
% containing a transcript identifier.
%
% This function is useful for exporting a transcriptome into a file format that
% the Bowtie read aligner software can use to build a new index.
%
% Inputs:
%     filename - Name of the file into which the data should be written.

% Author: Matti Annala <matti.annala@tut.fi>


function [] = write_seq_fasta(transcripts, filename)

output = fopen(filename, 'W');

if isfield(transcripts, 'name')
	for k = 1:length(transcripts.name)
		fprintf(output, '>%s\n%s\n', transcripts.name{k}, ...
			transcripts.sequence{k});
	end
elseif isfield(transcripts, 'sequence')
	for k = 1:length(transcripts.sequence)
		fprintf(output, '>%d\n%s\n', k, transcripts.sequence{k});
	end
elseif iscellstr(transcripts)
	for k = 1:length(transcripts)
		fprintf(output, '>%d\n%s\n', k, transcripts{k});
	end
else
	error 'Sequence structure is of unrecognized type.';
end

fclose(output);

