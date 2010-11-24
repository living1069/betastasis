function [genome, transcriptome, exons] = read_ncbi_refseq()

chromosomes = struct;
chromosomes.Name = {};
for k = 1:22, chromosomes.Name{k} = num2str(k); end
chromosomes.Name{23} = 'X';
chromosomes.Name{24} = 'Y';
chromosomes.Name{25} = 'M';

genes = struct;
genes.Name = {};
genes.Chromosome = [];
genes.TranscriptCount = [];
genes.Transcripts = [];

transcripts = struct;
transcripts.Name = {};
transcripts.Sequence = {};
transcripts.Gene = [];
transcripts.ExonCount = [];
transcripts.Exons = [];
transcripts.ExonPositions = [];

exons = struct;
exons.Transcript = [];
exons.Coordinates = [];

gene_count = 0;
transcript_count = 0;
exon_count = 0;

gene_map = containers.Map;
transcript_map = containers.Map;
exon_map = containers.Map;

assembly_len = 0;

fprintf(1, 'Building contig map...\n');
contig_map = build_contig_map()

fprintf(1, 'Parsing genomic annotations...\n');

gbff = fopen('refseq_human.genomic.gbff');
if gbff == -1
	fprintf(1, 'File %s could not be opened.', filepath);
	return;
end

parse_mode = 0;

while 1
	line = fgetl(gbff);
	if line == -1, break, end
		
	fprintf(1, '%s\n', line);
	
	while 1
		if parse_mode == 0
			if length(line) >= 7 && strcmp(line(1:7), 'PRIMARY')
				fprintf(1, '>> Parsing gene coordinates...\n');
				parse_mode = 1;
				assembly_len = 0;
				prev_frag_id = '';
				prev_frag_offset = 0;
				break;
			end
			
			matches = regexp(line, '^\s+gene\s+[<>]?(\d+)\.\.[<>]?(\d+)');
			if length(matches) == 1
				fprintf(1, '>> Parsing gene attributes...\n');
				parse_mode = 2;
				break;
			end
			
			break;
		
			
			
			
			
		elseif parse_mode == 1
			tokens = regexp(line, '^\s+\d+-\d+\s+(\S+)\s+(\d+)-(\d+)(.*)', ...
				'tokens');
			if length(tokens) ~= 1
				if length(line) < 8 || ~strcmp(line(1:8), 'FEATURES')
					line, error '>> Abrupt termination of PRIMARY coords.';
				end
				
				fprintf(1, '>> Stopped parsing gene coordinates...\n');
				parse_mode = 0;
				continue;
			end
			
			% The goal here is to parse the genome fragments and to flip
			% them so that all coordinates point to the + strand.
			
			tokens = tokens{1};
			frag_id = tokens{1};
			frag_start = str2num(tokens{2});
			frag_end = str2num(tokens{3});
			if strfind(tokens{4}, 'c')
				
			end
			
			if frag_id(1) == '"'
				% Either some bases have been changed or new bases have been
				% inserted.
				if assembly_len == 0
					error '>> Assembly starts with a hotfix. Not supported.'; 
				end
				
				fixed_seq = frag_id(2:end-1);
				prev_frag_offset = prev_frag_offset + length(fixed_seq);
				assembly_len = assembly_len + length(fixed_seq);
			else
				% The sequence is a component of some contig.
				if ~contig_map.isKey(frag_id)
					error('>> Contig %s was not found.', frag_id);
				end
				
				contig = contig_map(frag_id);
				if assembly_len == 0
					assembly_base = contig.Start - ...
						(frag_start - contig.CompStart) - 1;
					fprintf(1, '>> First fragment starts at %d.\n', ...
						assembly_base);
				elseif strcmp(frag_id, prev_frag_id)
					if frag_start ~= prev_frag_offset
						error '>> Gap in fragment.';
					end
					fprintf(1, '>> Previous fragment continues.\n');
				else
					if frag_start ~= contig.CompStart
						error '>> Invalid fragment start.';
					elseif prev_frag_offset ~= 0
						error '>> Previous fragment did not complete.';
					end
					fprintf(1, '>> Next fragment.\n');
				end
				
				if frag_end ~= contig.CompEnd
					prev_frag_offset = frag_end + 1;
				else
					prev_frag_offset = 0;
				end
				
				assembly_len = assembly_len + (frag_end - frag_start + 1);
				prev_frag_id = frag_id;
			end
			
			break;
		
			
			
			
			
		elseif parse_mode == 2
			tokens = regexp(line, '/gene="(.+)"', 'tokens');
			if length(tokens) == 1
				tokens = tokens{1}; gene_id = tokens{1};
				if gene_map.isKey(gene_id)
					error 'Gene defined twice.';
				end
				
				gene_count = gene_count + 1;
				gene_map(gene_id) = gene_count;
				genes.Name{gene_count} = gene_id;
				break;
			end
			
			if length(line) < 6 || (strcmp(line(1:5), '     ') && line(6) ~=' ')
				fprintf(1, '>> Stopped parsing gene attributes...\n');
				parse_mode = 0;
				continue;
			end
			
			break;

		end
	end
	
%	matches = regexp(line, '\s+exon\s+(\d+)\.\.(\d+)', 'tokens');
%	if length(matches) == 1
%		tokens = matches{1};
%		exon_transcripts(current_exon) = k;
%		exon_pos(current_exon, 1:2) = ...
%			[ str2double(tokens{1}), str2double(tokens{2}) ];
%		current_exon = current_exon + 1;
%		continue;
%	end
	
%	matches = regexp(line, '^\s+variation\s+[^\d]*(\d+)', 'tokens');
%	if length(matches) == 1
%		tokens = matches{1};
%		
%		snp_count = snp_count + 1;
%		snp_transcripts(snp_count) = k;
%		snp_pos(snp_count) = str2double(tokens{1});
%		
%		found_original = 0;
%		
%		while 1
%			snp_line = fgetl(rna_gbff);
%			if snp_line == -1, break, end
%			
%			matches = regexp(snp_line, '/replace="([tcga]*)"', 'tokens');
%			if length(matches) == 1
%				tokens = matches{1};
%				
%				if found_original
%					snp_variations{snp_count} = tokens{1};
%					%fprintf(1, 'Found SNP %s -> %s in transcript %d\n', ...
%					%	snp_originals{snp_count}, snp_variations{snp_count}, k);
%					break;
%				else
%					snp_originals{snp_count} = tokens{1};
%					found_original = 1;
%				end
%			elseif found_original
%				% Found a type of SNP that we don't recognize yet. Usually this
%				% means a "largedeletion" SNP.
%				%fprintf(1, 'Invalid SNP found in transcript %d.\n', k);
%				snp_count = snp_count - 1;
%				break;
%			end
%		end
%	end
end

fclose(gbff);
return;
	
	




function cmap = build_contig_map(chromosomes)
d = dir('.');
files = {d.name};

cmap = containers.Map;

for k = 1:length(files)
	file = files{k};
	tokens = regexpi(file, '.+chr(.+)\.agp', 'tokens');
	if length(tokens) ~= 1, continue, end
	
	chr = tokens{1}; chr = chr{1};
	chrnum = chromosome_sym2num(chr);
	
	fid = fopen(file);
	while 1
		line = fgetl(fid);
		if length(line) >= 6 && strcmp(line(1:6), '#chrom'), break, end
	end
	
	data = textscan(fid, '%*s %d %*d %*s %s %s %s %s %s', ...
		'Delimiter', '\t', 'ReturnOnError', 0);
	
	start = data{1};
	type = data{2};
	name = data{3};
	comp_start = data{4};
	comp_end = data{5};
	strand = data{6};
	clear data;
	
	for c = 1:length(name)
		if strcmp(type, 'N'), continue, end
		
		contig = struct;
		contig.Chromosome = chrnum;
		contig.Start = start(c);
		contig.CompStart = str2num(comp_start{c});
		contig.CompEnd = str2num(comp_end{c});
		contig.Strand = strand{c};
		
		cmap(name{c}) = contig;
	end
	
	fclose(fid);
end

return;




function [] = find_regions_and_hotfixes()

regions = struct;
regions.Name = {};
regions.Chromosome = [];
regions.Strand = true(10, 1);
regions.Start = [];
regions.End = [];

hotfixes = struct;
hotfixes.Chromosome = [];
hotfixes.Coordinate = [];

gbff = fopen('refseq_human.genomic.gbff');
if gbff == -1
	fprintf(1, 'File %s could not be opened.', filepath);
	return;
end

parse_mode = 0;

while 1
	line = fgetl(gbff);
	if line == -1, break, end
		
	fprintf(1, '%s\n', line);
	
	while 1
		if parse_mode == 0
			if length(line) >= 7 && strcmp(line(1:7), 'PRIMARY')
				fprintf(1, '>> Parsing region coordinates...\n');
				parse_mode = 1;
				assembly_len = 0;
				prev_frag_id = '';
				prev_frag_offset = 0;
				break;
			end
			
			matches = regexp(line, '^\s+gene\s+[<>]?(\d+)\.\.[<>]?(\d+)');
			if length(matches) == 1
				fprintf(1, '>> Parsing gene attributes...\n');
				parse_mode = 2;
			end
			
			break;
		
			
			
			
			
		elseif parse_mode == 1
			tokens = regexp(line, '^\s+\d+-\d+\s+(\S+)\s+(\d+)-(\d+)(.*)', ...
				'tokens');
			if length(tokens) ~= 1
				if length(line) < 8 || ~strcmp(line(1:8), 'FEATURES')
					line, error '>> Abrupt termination of PRIMARY coords.';
				end
				
				fprintf(1, '>> Stopped parsing gene coordinates...\n');
				parse_mode = 0;
				continue;
			end
			
			% The goal here is to parse the genome fragments and to flip
			% them so that all coordinates point to the + strand.
			
			tokens = tokens{1};
			frag_id = tokens{1};
			frag_start = str2num(tokens{2});
			frag_end = str2num(tokens{3});
			if strfind(tokens{4}, 'c')
				
			end
			
			if frag_id(1) == '"'
				% Either some bases have been changed or new bases have been
				% inserted.
				if assembly_len == 0
					error '>> Assembly starts with a hotfix. Not supported.'; 
				end
				
				fixed_seq = frag_id(2:end-1);
				prev_frag_offset = prev_frag_offset + length(fixed_seq);
				assembly_len = assembly_len + length(fixed_seq);
			else
				% The sequence is a component of some contig.
				if ~contig_map.isKey(frag_id)
					error('>> Contig %s was not found.', frag_id);
				end
				
				contig = contig_map(frag_id);
				if assembly_len == 0
					assembly_base = contig.Start - ...
						(frag_start - contig.CompStart) - 1;
					fprintf(1, '>> First fragment starts at %d.\n', ...
						assembly_base);
				elseif strcmp(frag_id, prev_frag_id)
					if frag_start ~= prev_frag_offset
						error '>> Gap in fragment.';
					end
					fprintf(1, '>> Previous fragment continues.\n');
				else
					if frag_start ~= contig.CompStart
						error '>> Invalid fragment start.';
					elseif prev_frag_offset ~= 0
						error '>> Previous fragment did not complete.';
					end
					fprintf(1, '>> Next fragment.\n');
				end
				
				if frag_end ~= contig.CompEnd
					prev_frag_offset = frag_end + 1;
				else
					prev_frag_offset = 0;
				end
				
				assembly_len = assembly_len + (frag_end - frag_start + 1);
				prev_frag_id = frag_id;
			end
			
			break;
			
		end
	end


return;






function [] = patch_genome(hotfixes)

return;
