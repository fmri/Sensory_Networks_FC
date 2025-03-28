% Reorder ROIs for conn_gPPI_res_aAvA.mat so that adjacency matrix is
% easier to create/use

load('ROI.mat');

ROI_data = cat(4, ROI.y); % subject x connection_ROI x condition x ROI
beta_diffs = squeeze(ROI_data(:,:,2,:) - ROI_data(:,:,1,:)); % subj x connection_ROI x ROI
beta_diffs = permute(beta_diffs, [1,3,2]); % subj x ROI x connection_ROI
ROI_names = ROI.names;

desired_order = {'tgPCS (L)', 'FO (L)', 'CO (L)', 'cIFSG (L)', 'cmSFG (L)', 'tgPCS (R)', 'FO (R)', 'CO (R)', ...
    'cIFSG (R)', 'cmSFG (R)', 'pAud (L)', 'pAud (R)',  ...
    'sPCS (L)', 'iPCS (L)', 'midIFS (L)', 'sPCS (R)', 'iPCS (R)', 'midIFS (R)', 'pVis (L)', 'pVis (R)', ...
    'aINS (L)', 'preSMA (L)', 'dACC (L)', 'aINS (R)', 'preSMA (R)', 'dACC (R)'};
ROI_str = 'avsm_ROIs';
names_clean = cellfun(@(f) f(length(ROI_str)+2:end), ROI_names, 'UniformOutput', false);
[ROIs_match, reorder_inds] = ismember(desired_order, names_clean);

ROI_names = names_clean(reorder_inds);
beta_diffs = beta_diffs(:,reorder_inds, reorder_inds);

% Make adjacency matrix 
adj_mat = zeros(26);
adj_mat(1:10, 1:10) = 1; % frontal auditory
adj_mat(11:12, 11:12) = 1; % p aud
adj_mat(13:18, 13:18) = 1; % frontal visual
adj_mat(19:20, 19:20) = 1; % p vis
adj_mat(21:26, 21:26) = 1; % MD


%save('conn_gPPI_betadiffs.mat', "beta_diffs", "ROI_names", "adj_mat");
save('conn_gPPI_betadiffs.mat', "beta_diffs", "ROI_names");



