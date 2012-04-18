classdef BowtieAlignments
	properties (Hidden = true, Access = private)
		fid = []
		batch = 1e6
		fields = {}
		prev_idx = 0
		num_alignments = 0
	end
	
	methods
		function obj = BowtieAlignments(url, num_al, varargin)
			
			for k = 1:2:length(varargin)
				if strcmpi(varargin{k}, 'Fields')
					obj.fields = varargin{k+1};
					continue;
				end
				
				if regexpi(varargin{k}, 'Batch')
					obj.batch = varargin{k+1};
					continue;
				end
				
				error('Unrecognized option "%s".', varargin{k});
			end
		
			obj.fid = fopen(url);
			obj.num_alignments = num_al;
		end
		
		function [m, n] = size(obj)
			m = 1;
			n = ceil(obj.num_alignments / obj.batch);
			if obj.batch == Inf, n = 1; end
		end
		
		function al = subsref(obj, s)
			if length(s) ~= 1 || ~strcmp(s(1).type, '()')
				error 'Invalid access.';
			end
			
			format = '%s%s%s%d%s';
			if obj.batch == Inf
				data = textscan(obj.fid, format, ...
					'Delimiter', '\t', 'ReturnOnError', false);
			else
				data = textscan(obj.fid, format, obj.batch, ...
					'Delimiter', '\t', 'ReturnOnError', false);
			end
			
			al = struct;
			al.read = data{1};
			al.strand = char(data{2});
			al.target = data{3};
			al.offset = data{4};
			al.sequence = data{5};
		end
	end
end
