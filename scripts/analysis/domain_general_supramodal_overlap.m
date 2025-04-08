%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to assess the overlap between supramodal ROIs,
% and Assem et al. 2022 domain general ROIs using dice coefficients.
% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/'));
ccc;

%% Initialize Key Variables
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'}; % 21 subjs - union of resting state and localizer analyses
N = length(subjCodes);

ROI_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs/';
sm_ROI_names = {'sm_sPCS', 'sm_iPCS', 'sm_midFSG', 'sm_aINS', 'sm_dACC', 'sm_preSMA'};
N_sm = length(sm_ROI_names);

dg_ROI_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/surf/';
DG_ROI_names = {'s6-8', 'i6-8', '8C', 'IFJp', 'p9-46v', '6r', 'p10p', 'a10p', '11l',...
        'a47r', 'p47r', 'a9-46v', 'FOP5', 'AVI', 'TE1m', 'TE1p', 'AIP', 'IP2',...
        'IP1', 'LIPd', 'MIP', 'PGs', 'PFm', 'POS2', 'SCEF', '8BM', 'a32pr', 'd32'};
N_dg = length(DG_ROI_names);

hemis = {'lh', 'rh'};   
N_hemis = length(hemis);

overlap = nan(N, N_hemis, N_dg, N_sm);
dice_coeffs = nan(N, N_hemis, N_dg, N_sm);
avsm_poverlap = nan(N, N_hemis, N_dg, N_sm);
dg_poverlap = nan(N, N_hemis, N_dg, N_sm);

%% Loop through subjs and ROIs and assess overlap 
for hh = 1:N_hemis
    hemi = hemis{hh};

    % Load DG ROIS
    [~, labels_dg, ctable_dg] = read_annotation([dg_ROI_dir hemi '.HCP-MMP1.annot']); 

    for nn = 1:N
        subjCode = subjCodes{nn};
    
        % Load annotation files
        [~, labels_sm, ctable_sm] = read_annotation([ROI_dir hemi '.' subjCode '_avsm_ROIs.annot']); 

        for sm = 1:N_sm
            sm_ROI_name = sm_ROI_names{sm};
            sm_ROI_code = ctable_sm.table(ismember(ctable_sm.struct_names, sm_ROI_name),5);
            sm_ROI_mask = labels_sm == sm_ROI_code;
            for dg = 1:N_dg
                dg_ROI_name = [upper(hemi(1)) '_' DG_ROI_names{dg} '_ROI'];
                dg_ROI_code = ctable_dg.table(ismember(ctable_dg.struct_names, dg_ROI_name),5);
                dg_ROI_mask = labels_dg == dg_ROI_code;
                
                overlap(nn,hh,dg,sm) = sum(dg_ROI_mask & sm_ROI_mask);
                dice_coeffs(nn,hh,dg,sm) = (2*sum(dg_ROI_mask & sm_ROI_mask)) / (sum(dg_ROI_mask) + sum(sm_ROI_mask));
                avsm_poverlap(nn,hh,dg,sm) = sum(dg_ROI_mask & sm_ROI_mask)/sum(sm_ROI_mask);
                dg_poverlap(nn,hh,dg,sm) = sum(dg_ROI_mask & sm_ROI_mask)/sum(dg_ROI_mask);
            end

        end
    end
end


%% Plot heatmap of dice coef matrix
hem_avg_dice = squeeze(mean(dice_coeffs, 2))
groupavg_dice= squeeze(mean(dice_coeffs, [1,2]));

groupavg_overlap= squeeze(mean(avsm_poverlap, [1,2]));

groupavg_dg_overlap= squeeze(mean(dg_poverlap, [1,2]));


figure; 
heatmap(replace(sm_ROI_names, '_', ' '), DG_ROI_names, groupavg_overlap)
title('dice');

figure; 
heatmap(replace(sm_ROI_names, '_', ' '), DG_ROI_names, groupavg_overlap)
title('sm overlap');

figure; 
heatmap(replace(sm_ROI_names, '_', ' '), DG_ROI_names, groupavg_dg_overlap)
title('dg overlap');

%% Test Dice significant difference from 0
ps = nan(N_dg, N_sm);
ts = nan(N_dg, N_sm);
sds = nan(N_dg, N_sm);
for dg = 1:N_dg
    for sm = 1:N_sm
        [~,ps(dg,sm),~,t] = ttest(squeeze(hem_avg_dice(:,dg,sm)));
        ts(dg,sm) = t.tstat;
        sds(dg,sm) = t.sd;
    end
end

group_sd_nohem = squeeze(std(hem_avg_dice,1));

figure; 
heatmap(replace(sm_ROI_names, '_', ' '), DG_ROI_names, group_sd_nohem)
title('Dice SDs');

figure; 
heatmap(replace(sm_ROI_names, '_', ' '), DG_ROI_names, ps);
title('Dice pvals');
