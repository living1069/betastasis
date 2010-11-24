function gene_expr = gene_expression_from_transcript_expression( ...
	transcript_expr, summarization_method)

global organism;
genome = organism.Genes;

gene_expr = struct('Mean', zeros(length(genome.Name), ...
	size(transcript_expr.Mean, 2)));
	
for g = 1:length(genome.Name)
	gene_transcripts = genome.Transcripts(g, 1:genome.TranscriptCount(g));
	gene_transcripts = gene_transcripts( ...
		~isnan(sum(transcript_expr.Mean(gene_transcripts, :), 2)));
	
	if length(gene_transcripts) == 0
		gene_expr.Mean(g, :) = NaN;
		continue;
	end
	
	if strcmpi(summarization_method, 'sum')
		gene_expr.Mean(g, :) = ...
			sum(transcript_expr.Mean(gene_transcripts, :), 1);
	elseif strcmpi(summarization_method, 'mean')
		gene_expr.Mean(g, :) = ...
			mean(transcript_expr.Mean(gene_transcripts, :), 1);
	elseif strcmpi(summarization_method, 'max')
		gene_expr.Mean(g, :) = ...
			max(transcript_expr.Mean(gene_transcripts, :), [], 1);
	else
		error 'Unrecognized summarization method requested.';
	end
end

