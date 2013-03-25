function genes = transcript_to_hugo(transcripts, refgene_file)

data = readtable(refgene_file, 'Header', false, 'Include', [2 13]);
ucsc_tx = data{1};
ucsc_hugo = data{2};

[valid, pos] = ismember(transcripts, ucsc_tx);
genes = repmat({''}, length(transcripts), 1);
genes(valid) = ucsc_hugo(pos(valid));

