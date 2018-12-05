clear all;

%% read bul_rs 
% n = 100, w/ motion thresh mean motion<0.5, max motion<4.5

load('all_mats_rs.mat')
i = [1	2	3	4	5	6 7	8	9	10		11		12	13	14		15	16	17	18	19		20	21		22	23	24		25	26	27	28	29	30			31	32	33			34	35	36	37	38	39	40	41	42	43	44		45	46	47	48		49	50	51	52	53	54					55	56		57		58	59	60		61	62	63	64	65		66	67	68	69	70	71		72	73		74		75	76		77	78	79	80	81	82	83	84	85	86	87	88	89	90		91	92	93	94	95	96		97]
all_mats_rs = allmat(:,:,i);
    
%% load selected overlapping edges from predictive model
% load ('100_overlap.mat') 
% %load ('90_overlap.mat')

no_node = 227 %shen, cereb removed
no_subs = 97

pos_mask = zeros(no_node, no_node);
neg_mask = zeros(no_node, no_node);

pos_edges = overlap_mask_pos;
neg_edges = overlap_mask_neg;

pos_mask(pos_edges) = 1;
neg_mask(neg_edges) = 1; 

% soc mask
soc_mask = xlsread("E:\!Study\study_cpm_pols\external_vals\hk_cpm_pols_replication\ormask.csv");
pos_mask = pos_mask .* soc_mask;
neg_mask = neg_mask .* soc_mask;

%% sum all edges using pos/neg mask from ts FC
all_mats_rs(isnan(all_mats_rs))=0;

train_mats = all_mats_rs;

train_sumpos = zeros(97,1); %train mat size external validation sample size·Î ¼öÁ¤
train_sumneg = zeros(97,1); 

    for ss = 1:size (train_sumpos)
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask))/2;
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask))/2;
    end

