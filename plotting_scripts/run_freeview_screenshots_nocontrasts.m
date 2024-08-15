%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to call freeview_screenshots to take
% screenshots of all ROIs on all subjs and save them.
%
% Tom Possidente August 1 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'))
ccc;

%% Get subjIDs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'AH', 'SL'})); % no ROIs for these subjs
N = length(subjCodes);

use_fsaverage = true; % if false will use the ROIs from individual surface

%% Create lookup table for ROIs and corresponding colors
ROI = ["aINS", "preSMA", "ppreCun", "dACC", ... % multisensory
        "sPCS", "iPCS", "midIFS", "aIPS", "pIPS", "DO", "LOT", "VOT",... % visual
        "tgPCS", "cIFSG", "pAud", "CO", "FO", "cmSFG"]'; % auditory
% color = [[0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; ... % blue
%           [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; ... % green
%           [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]]; % green
color = [[0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; ... % green
          [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; ... % blue
          [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]]; % orange
ctable = table(ROI, color(:,1), color(:,2), color(:,3), 'VariableNames', {'ROI', 'c1', 'c2', 'c3'});


%% Loop through subjs and get ROI paths
lh_ROI_paths = {};
rh_ROI_paths = {};

roi_base_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';

all_roi_files = {dir(roi_base_dir).name};

if use_fsaverage
    all_roi_files = all_roi_files(~contains(all_roi_files, 'indiv'));
    ss_suffix = '_fsavg';
else
    all_roi_files = all_roi_files(contains(all_roi_files, 'indiv'));
    ss_suffix = '_self';
end


for ss = 1:N % loop through subjs
    subjCode = subjCodes{ss};
    subj_ROIs = all_roi_files(contains(all_roi_files, [subjCode '_'])); % find all ROI files for subj
    if strcmp(subjCode, 'NS') % due to "NS" being in the label aINS, we have to take out the wrong aINS for this subj
        bad_aINS = contains(subj_ROIs, 'aINS') & ~contains(subj_ROIs, 'NS_aINS');
        subj_ROIs = subj_ROIs(~bad_aINS);
    end

    lh_ROI_files = subj_ROIs(contains(subj_ROIs, 'lh')); % find all lh ROIs
    rh_ROI_files = subj_ROIs(contains(subj_ROIs, 'rh')); % find all rh ROIs
    

    for rr = 1:length(lh_ROI_files) % for each ROI
        ROI_split = split(lh_ROI_files{rr}, '_');
        ROI_name = ROI_split{2}; % isolate name of ROI only 
        lh_ROI_paths{ss, 1}{rr} = [roi_base_dir lh_ROI_files{rr}]; % put it in correct cell for subj and contrast
    end

    for rr = 1:length(rh_ROI_files)
        ROI_split = split(rh_ROI_files{rr}, '_');
        ROI_name = ROI_split{2};
        rh_ROI_paths{ss, 1}{rr} = [roi_base_dir rh_ROI_files{rr}];
    end

end

%% Actually run freeview_screenshots()
contrasts = {};
freeview_screenshots(subjCodes, contrasts, lh_ROI_paths, rh_ROI_paths, ctable, use_fsaverage,[],[],[],[],[],[],[],'0.9',ss_suffix)
crop_ppt_fv_images(subjCodes, contrasts, true, [], 'ROI_fsavg_surface_screenshots_nocontrasts')
