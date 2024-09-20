%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to call freeview_screenshots to take
% screenshots of all ROIs on the space time task contrasts for
% all subjs and save them.
%
% Tom Possidente September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'))
ccc;

%% Get subjIDs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'AH', 'SL', 'RR'})); % no ROIs for these subjs
N = length(subjCodes);

%% Create lookup table for ROIs and corresponding colors
ROIs = ["aINS", "preSMA", "ppreCun", "dACC", ... % multisensory
    "sPCS", "iPCS", "midIFS", "pVis_mod",... % visual
    "tgPCS", "cIFSG", "pAud", "CO", "FO", "cmSFG"]'; % auditory
N_ROIs = length(ROIs);
color = [[0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; ... % blue
    [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; ... % green
    [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]]; % green
% color = [[0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; ... % green
%           [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; ... % blue
%           [256 165 0]; [256 165 0]; [256 165 0]; [256 165 0]; [256 165 0]; [256 165 0]]; % orange
ctable = table(ROIs, color(:,1), color(:,2), color(:,3), 'VariableNames', {'ROI', 'c1', 'c2', 'c3'});


%% Set ROI and contrast lists
%contrast_list = {'A-P', 'V-A'};
AV_A_P_ROIs = ["aINS", "preSMA", "ppreCun", "dACC"];
sV_pV_ROIs = ["sPCS", "iPCS", "midIFS", "pVis_mod"];
tA_pA_ROIs = ["tgPCS", "cIFSG", "pAud", "CO", "FO", "cmSFG"];

lh_ROI_paths = cell(N, 3);
rh_ROI_paths = cell(N, 3);
contrasts = cell(N,3);

roi_base_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';

all_roi_files = {dir(roi_base_dir).name};

for ss = 1:N % loop through subjs
    subjCode = subjCodes{ss};
    subj_ROIs = all_roi_files(contains(all_roi_files, subjCode) & ~contains(all_roi_files,'indiv')); % find all ROI files for subjs
    if strcmp(subjCode, 'NS') % due to "NS" being in the label "aINS", we have to take out the wrong aINS for this subj
        bad_aINS = contains(subj_ROIs, 'aINS') & ~contains(subj_ROIs, 'NS_aINS');
        subj_ROIs = subj_ROIs(~bad_aINS);
    end
    lh_ROI_files = subj_ROIs(contains(subj_ROIs, 'lh')); % find all lh ROIs
    rh_ROI_files = subj_ROIs(contains(subj_ROIs, 'rh')); % find all rh ROIs

    contrasts{ss,1} = 'AV_A-P';
    contrasts{ss,2} = 'sV-pV';
    contrasts{ss,3} = 'tA-pA';

    for rr = 1:N_ROIs % for each ROI
        fname_suffix = '';
        ROI_curr = lh_ROI_files(contains(lh_ROI_files, ROIs{rr}) & contains(lh_ROI_files, '.label'));
        if isempty(ROI_curr)
            disp(['Subj ' subjCode ' missing ROI ' ROIs{rr} ' lh']);
            continue;
        end
        if length(ROI_curr) > 1 || contains(ROI_curr, 'replacement')
            ROI_curr = ROI_curr(contains(ROI_curr, 'replacement'));
            assert(~isempty(ROI_curr));
            fname_suffix = '_replacement';
        end

        ROI_split = split(ROI_curr, '_');
        ROI_name = ROI_split{2}; % isolate name of ROI only
        if ismember(ROI_name, AV_A_P_ROIs) % see which contrast this ROI belongs to
            cc = 1;
        elseif ismember(ROI_name, sV_pV_ROIs)
            cc = 2;
        elseif ismember(ROI_name, tA_pA_ROIs)
            cc = 3;
        else
            disp(['Subj ' subjCode ' ROI not used: ' ROI_name]);
            continue;
        end
        lh_ROI_paths{ss, cc}{end+1} = [roi_base_dir subjCode '_' ROI_name '_lh' fname_suffix '.label']; % put it in correct cell for subj and contrast
        assert(~isempty(lh_ROI_paths{ss, cc}{end}))
    end

    for rr = 1:N_ROIs % for each ROI
        fname_suffix = '';
        ROI_curr = rh_ROI_files(contains(rh_ROI_files, ROIs{rr}) & contains(rh_ROI_files, '.label'));
        if isempty(ROI_curr)
            disp(['Subj ' subjCode ' missing ROI ' ROIs{rr} ' rh']);
            continue;
        end
        if length(ROI_curr) > 1 || contains(ROI_curr, 'replacement')
            ROI_curr = ROI_curr(contains(ROI_curr, 'replacement') );
            assert(~isempty(ROI_curr));
            fname_suffix = '_replacement';
        end
        ROI_split = split(ROI_curr, '_');
        ROI_name = ROI_split{2}; % isolate name of ROI only
        if ismember(ROI_name, AV_A_P_ROIs) % see which contrast this ROI belongs to
            cc = 1;
        elseif ismember(ROI_name, sV_pV_ROIs)
            cc = 2;
        elseif ismember(ROI_name, tA_pA_ROIs)
            cc = 3;
        else
            disp(['Subj ' subjCode ' ROI not used: ' ROI_name]);
            continue;
        end
        rh_ROI_paths{ss, cc}{end+1} = [roi_base_dir subjCode '_' ROI_name '_rh' fname_suffix '.label']; % put it in correct cell for subj and contrast
        assert(~isempty(rh_ROI_paths{ss, cc}{end}))
    end

end

%% Actually run freeview_screenshots()
func_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/';
func_folder = 'bold';
analysis_name = 'spacetime_contrasts';
saveDir = '/projectnb/somerslab/tom/projects/spacetime_network/figures_images/spacetime_contrasts/';
freeview_screenshots(subjCodes, contrasts, lh_ROI_paths, rh_ROI_paths, ctable, [], saveDir, func_path, func_folder, [], ...
    analysis_name, [], [], '0.5');
crop_ppt_fv_images(subjCodes, contrasts, [], saveDir, 'spacetime_contrasts.pptx')
