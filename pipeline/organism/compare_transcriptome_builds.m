
% This function performs a comparison on two different transcriptome builds.
% The comparison is based on transcript identifiers, which means that the
% compared transcriptomes *must* use the same transcript naming convention.
% If the naming conventions differ, all transcripts are shown as having changed.
%
% Inputs:
%     ref_transcriptome - Reference transcriptome build. All reported
%         transcriptome changed are relative to this build.
%     new_transcriptome - New transcriptome build.
%
% Outputs:
%     transcriptome_changes - Structure that contains a list of added, removed
%         and changed transcripts relative to the reference build.
%
% Author: Matti Annala <matti.annala@tut.fi>

function transcriptome_changes = compare_transcriptome_builds(...
	ref_transcriptome, new_transcriptome)
	
added = {};
removed = {};
changed = {};

old_keys = keys(old_genes);

for k = 1:length(old_keys)
	if isKey(new_genes, old_keys(k))
		sequence = new_genes(old_keys{k});
		if ~strcmp(sequence, old_genes(old_keys{k}))
			changed{length(changed) + 1} = old_keys{k};
		end
	else
		removed{length(removed) + 1} = old_keys{k};
	end
end

new_keys = keys(new_genes);

for k = 1:length(new_keys)
	if ~isKey(old_genes, new_keys{k})
		added{length(added) + 1} = new_keys{k};
	end
end

fprintf(1, 'Transcriptome comparison results:\n');
fprintf(1, '- %d transcripts added\n', length(added));
fprintf(1, '- %d transcripts removed\n', length(removed));
fprintf(1, '- %d transcripts changed\n', length(changed));

transcriptome_changes = struct( ...
	'Added', added, ...
	'Removed', removed, ...
	'Changed', changed);

