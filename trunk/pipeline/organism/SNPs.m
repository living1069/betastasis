classdef SNPs < handle
	properties (Hidden = true, Access = private)
		Chromosomes
	end
	
	methods
		function obj = SNPs(orgname, version)
			global organism;
			obj.Chromosomes = cell(length(organism.Chromosomes.Name), 1);
		end
		
		function ret = subsref(obj, s)
			global organism;
			
			if ~strcmp(s(1).type, '{}'), error 'Invalid access.'; end
			
			idx = s(1).subs;
			if length(idx) ~= 1, error 'Invalid access.'; end
			chr = idx{1};
			
			if isempty(obj.Chromosomes{chr})
				snps = load([ppath '/organisms/' flatten_str(organism.Name) ...
					'/' flatten_str(organism.Version) ...
					'/snps/chr' organism.Chromosomes.Name{chr} '.mat']);
				obj.Chromosomes{chr} = snps.snps;
			end
			
			ret = obj.Chromosomes{chr};
			
			if length(s) > 1
				ret = subsref(ret, s(2:end));
			end
		end
	end
end


