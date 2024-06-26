%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to make a heatmap plot of the resting state
% functional connectivity between ROIs for all available subjs using Conn's
% precomputed GLM connectivity analysis results
% Tom Possidente - June 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccc;

%%
rs_conn_results_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_resting_state/results/firstlevel/RRC_01/resultsROI_Condition001.mat';
load(rs_conn_results_path, 'Z','names');
avg_connmat = mean(Z(:,:,[1,2,3,4,5,6,7,8,9,10,11,12,13]), 3);

aud_ROIs_use = {'tgPCS', 'cIFS/G', 'FO', 'CO', 'pAud'};
vis_ROIs_use = {'pVis', 'preSMA-V', 'SPCS', 'IPCS', 'midIFS'};
mult_ROIs_use = {'cmSFG_mult', 'Ins_mult'};

names_clean = cellfun(@(f) f(10:end), names, 'UniformOutput', false);
desired_order = {'tgPCS (L)', 'FO (L)', 'CO (L)', 'cIFS_G (L)', 'pAud (L)', 'tgPCS (R)', 'FO (R)', 'CO (R)', ...
                 'cIFS_G (R)', 'pAud (R)',...
                 'SPCS (L)', 'IPCS (L)', 'midIFS (L)', 'pVis (L)', 'SPCS (R)', 'IPCS (R)', 'midIFS (R)', 'pVis (R)',...
                 'Ins_mult (R)'};

[~, name_inds] = ismember(desired_order, names_clean);

connmat_reorder = avg_connmat(name_inds,name_inds);

%% Plot
%[min(Z(:,:,pp),[], 'all'),max(Z(:,:,pp),[],'all')]
for pp = 1:size(Z,3)
    figure;
    heatmap(desired_order, desired_order, Z(:,:,pp), 'Colormap', turbo, 'ColorLimits', [-1,1])
    title(['Conn Mat Subj ', num2str(pp)])
end

figure;
%[min(connmat_reorder,[], 'all'),max(connmat_reorder,[],'all')]
heatmap(desired_order, desired_order, connmat_reorder, 'Colormap', turbo, 'ColorLimits', [-1,1])
title('Group Average Conn Mat')