classdef Sequence
	properties (Hidden = true, Access = private)
		Filename
	end
	
	methods
		function obj = Sequence(seq_file)
			obj.Filename = seq_file;
		end
		
		function ret = subsref(obj, s)
			global organism;
			
			fid = fopen([ppath '/organisms/' flatten_str(organism.Name) ...
				'/' flatten_str(organism.Version) ...
				'/genome/' obj.Filename]);
			if numel(s) > 1 || ~strcmp(s.type, '()') || numel(s.subs) > 1
				error 'Invalid sequence access.';
			end
			
			idx = s.subs{1};
			idx_min = min(idx);
			idx_max = max(idx);
			
			fseek(fid, idx_min - 1, 'bof'); 
			full_seq = fread(fid, idx_max - idx_min + 1, 'char');
			fclose(fid);
			
			ret = char(full_seq(idx - idx_min + 1))';
		end
	end
end


