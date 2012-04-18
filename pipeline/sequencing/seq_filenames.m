
% Author: Matti Annala <matti.annala@tut.fi>

function filenames = seq_filenames(reads)

filenames = cell(length(reads), 1);

for s = 1:length(reads.url)
	if rx(reads.format{s}, 'FASTA')
		filenames{s} = {[reads.url{s} '.fa']};
	elseif rx(reads.format{s}, 'FASTQ')
		filenames{s} = {[reads.url{s} '.fq']};
	else
		error 'Unrecognized read format.';
	end
		
	if rx(reads.paired{s}, 'Paired')
		filenames{s} = cat(2, ...
			regexprep(filenames{s}, '(.+)(\..+?)$','$1_1$2'), ...
			regexprep(filenames{s}, '(.+)(\..+?)$', '$1_2$2'));
	end
	
	if rx(reads.format{s}, 'gzip')
		filenames{s} = strcat(filenames{s}, '.gz');
	end
end


