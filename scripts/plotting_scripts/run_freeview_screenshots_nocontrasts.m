%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to call freeview_screenshots to take
% screenshots of all ROIs on all subjs and save them.
%
% Tom Possidente August 1 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/'))
ccc;

%% Get subjIDs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'AH', 'SL', 'RR'})); % no ROIs for these subjs
N = length(subjCodes);

use_fsaverage = true; % if false will use the ROIs from individual surface

%% Create lookup table for ROIs and corresponding colors
ROI = ["sm_aINS", "sm_preSMA", "sm_dACC", "sm_sPCS", "sm_iPCS", "sm_midFSG" ... % supramodal
    "sPCS", "iPCS", "midIFS", "pVis",... % visual
    "tgPCS", "cIFSG", "pAud", "CO", "FO", "cmSFG"]'; % auditory
% color = [[0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; ... % blue
%           [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; ... % green
%           [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]]; % green
color = [[0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; ... % green
    [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; ... % blue
    [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]]; % orange
ctable = table(ROI, color(:,1), color(:,2), color(:,3), 'VariableNames', {'ROI', 'c1', 'c2', 'c3'});


%% Loop through subjs and get ROI paths
lh_ROI_paths = {};
rh_ROI_paths = {};

roi_base_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs_final/';

if use_fsaverage
    roi_base_dir = [roi_base_dir 'fsaverage/'];
    ss_suffix = '_fsavg';
else
    roi_base_dir = [roi_base_dir 'subjspace/'];
    ss_suffix = '_self';
end

for ss = 1:N % loop through subjs
    subjCode = subjCodes{ss};
    subj_dir = [roi_base_dir subjCode '/'];
    subj_ROIs = {dir(subj_dir).name};

    lh_ROI_files = subj_ROIs(contains(subj_ROIs, 'lh')); % find all lh ROIs
    rh_ROI_files = subj_ROIs(contains(subj_ROIs, 'rh')); % find all rh ROIs
    
    count = 0;
    for rr = 1:length(lh_ROI_files) % for each ROI
        ROI_name = split(lh_ROI_files{rr}, '.');
        if ismember(ROI_name{2}, ROI)
            count = count + 1;
            lh_ROI_paths{ss, 1}{count} = [subj_dir lh_ROI_files{rr}]; % put it in correct cell for subj and contrast
        end
    end

    count = 0;
    for rr = 1:length(rh_ROI_files)
        ROI_name = split(rh_ROI_files{rr}, '.');
        if ismember(ROI_name{2}, ROI)
            count = count + 1;
            rh_ROI_paths{ss, 1}{count} = [subj_dir rh_ROI_files{rr}]; % put it in correct cell for subj and contrast
        end    
    end

end

%% Actually run freeview_screenshots()
contrasts = {};
freeview_screenshots(subjCodes, contrasts, lh_ROI_paths, rh_ROI_paths, ctable, use_fsaverage,[],[],[],[],[],[],[],'0.9',ss_suffix)
crop_ppt_fv_images(subjCodes, contrasts, true, [], 'ROI_fsavg_surface_screenshots_nocontrasts')


%    subjIDs, contrast_list, lh_label_list, rh_label_list, colortable, use_fsaverage, ...
%    save_dir, func_path, func_folder, recon_dir, analysis_name, stat_file, overlayThreshold, label_opacity, ss_suffix, ss_on, overlay_path_override
