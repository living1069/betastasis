function meta = tampere_pca_clinical(meta)

meta.sample_id = map_sample_ids(meta.sample_id);

% First we check if the sample codes are in the latest format.
% If not, we convert them.
%meta = meta_tabular( ...
%	'/data/csb/datasets/tampere_pca/tampere_pca_clinical.txt', meta);





function sample_id = map_sample_ids(sample_id)

regexps = { ...
 '8131', 'PC_8131'; ...
 '226', 'PC_22603'; ...
 '4906', 'PC_4906'; ...
 '5934', 'PC_5934'; ...
 '9324', 'PC_9324'; ...
 '10286', 'PC_10286'; ...
 '24173', 'PC_24173'; ...
 '4980|PPC_150311_4', 'PC_4980'; ...
 '4538', 'PC_4538'; ...
 '6102', 'PC_6102'; ...
 '12517', 'PC_12517'; ...
 '13943', 'PC_13943'; ...
 '18307', 'PC_18307'; ...
 '20873', 'PC_20873'; ...
 '15194|PPC_150311_3', 'PC_15194'; ...
 '194', 'PC_19403'; ...
 '4786', 'PC_4786'; ...
 '6342', 'PC_6342'; ...
 '6488', 'PC_6488'; ...
 '15420', 'PC_15420'; ...
 '15760', 'PC_15760'; ...
 '17163', 'PC_17163'; ...
 '22392', 'PC_22392'; ...
 '470', 'PC_470'; ...
 '6174', 'PC_6174'; ...
 '7875', 'PC_7875'; ...
 '8438', 'PC_8438'; ...
 '14670', 'PC_14670'; ...
 '6864', 'PC_6864'; ...
 '17447', 'PC_17447'; ...
 '278[_.]2', 'CRPC_278.2'; ...
 '278', 'CRPC_278'; ...
 '305', 'CRPC_305'; ...
 '348', 'CRPC_348'; ...
 '489[_.]2', 'CRPC_489.2'; ...
 '489', 'CRPC_489'; ...
 '530', 'CRPC_530'; ...
 '531', 'CRPC_531'; ...
 '539', 'CRPC_539'; ...
 '541', 'CRPC_541'; ...
 '542', 'CRPC_542'; ...
 '543[_.]2', 'CRPC_543.2'; ...
 '543', 'CRPC_543'; ...
 '697', 'CRPC_697'; ...
 '261|TURP_150311_7', 'CRPC_261'; ...
 '435|TURP_150311_8', 'CRPC_435'; ...
 '651', 'BPH_651'; ...
 '659', 'BPH_659'; ...
 '665', 'BPH_665'; ...
 '677', 'BPH_677'; ...
 '689', 'BPH_689'; ...
 '701', 'BPH_701'; ...
 '688|PPC_150311_2', 'BPH_688'; ...
 '656|TURP_150311_5', 'BPH_656'; ...
 '671|TURP_150311_6', 'BPH_671'; ...
 '337', 'BPH_337'; ...
 '456', 'BPH_456'; ...
 '652|PPC_150311_1', 'BPH_652' };

for s = 1:length(sample_id)
	id = sample_id{s};
	for r = 1:size(regexps, 1)
		if rx(id, regexps{r, 1})
			sample_id{s} = regexps{r, 2}; break;
		end
	end
end
			
