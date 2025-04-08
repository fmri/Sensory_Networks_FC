%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to plot a measure of deviation from the
% group average connectivity matrix against average motion to see if there
% is a coarse relationship between movement and individual variability in connectivity 
% Tom Possidente - June 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccc;


%% Load in average motion data and connectivity data
motion_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/conn_toolbox_folder/conn_resting_state/results/qa/QA_2024_06_21_093046259/';
conndata_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/';
load([conndata_dir 'rs_group_connmat.mat'], 'names', 'connmat_group_nandiag');
load([conndata_dir 'rs_indiv_connmats.mat'], 'connmats_nandiag');

N = 13;

avg_motion = NaN(N,1);
connmat_corrs = NaN(N,1);

for ss = 1:N
    
    subjnum = sprintf( '%03d', ss ) ;
    load([motion_dir 'QA_COV.subject' subjnum '.mat'], 'results_info');
    avg_motion(ss) = results_info.Values(4);
    
    ind_connmat = connmats_nandiag(:,:,ss);
    connmat_corrs(ss) = corr(ind_connmat(:), connmat_group_nandiag(:), 'rows', 'complete');

end

%% Plot
figure;
scatter(avg_motion, connmat_corrs, 'filled');
xlabel('Average Motion (mm)');
ylabel('Correlation between Ind and Group Connectivity Matrices');
title('Resting State Motion vs. Connectivity QC');