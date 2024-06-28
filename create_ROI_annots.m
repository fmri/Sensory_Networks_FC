%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to create a valid annotation
% file with all ROIs for each subj, which Conn can then use
% Tom Possidente - June 2024
%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Set up directories and subj info
experiment_name = 'spacetime';
projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

% Get subj codes
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;

% Set ROIs to use and get ROI file list
ROIs_not_used = {'pSTS_lh', 'cmSFG_lh', 'S2_lh', 'midIns_lh', 'S1_lh', 'pSTS-Tac_lh', 'ppTac_lh', 'PMv_lh', 'posCingSulc_lh',...
    'IPCS_mult_lh', 'Lat_par_mult_lh', 'LatPar_mult_lh', 'SPCS_mult_lh', 'midIFS_mult_lh', 'pSTS_rh', 'cmSFG_rh',...
    'S2_rh', 'midIns_rh', 'S1_rh', 'pSTS-Tac_rh', 'ppTac_rh', 'PMv_rh', 'posCingSulc_rh', 'IPCS_mult_rh', ...
    'Lat_par_mult_rh', 'LatPar_mult_rh', 'SPCS_mult_rh', 'midIFS_mult_rh'};
all_ROIs_use = {'tgPCS', 'cIFS_G', 'FO', 'CO', 'pAud', 'pVis', 'preSMA-V', 'SPCS', ...
                'IPCS', 'midIFS', 'cmSFG_mult', 'Ins_mult'};

ROI_file_list_pre = {dir([projectDir '/data/ROIs/']).name};
ROI_file_list = ROI_file_list_pre(contains(ROI_file_list_pre, 'nooverlap')); % only taking ROI files that are binarized and without overlap

fs_number = 163842;
missing_ROIs_allsubj = cell(length(subjCodes),1);

% Read in reference annot file so we can keep the same structure for our annot files
annotpath = [projectDir 'data/recons/fsaverage/label/lh.aparc.annot'];
[verts_ref, labels_ref, ctable_ref] = read_annotation(annotpath); % Read in lh.aparc.annot file for reference annotation file structure

% Create vertices for all annot files to use (same number of vertices for all bc fsaverage space
annot_verts = verts_ref;

%% Start looping through subjs
for ss = 1:length(subjCodes)

    % Find all ROIs for subj
    subjCode = subjCodes{ss};
    file_mask = contains(ROI_file_list, [lower(subjCode) '_ROI']) & ~contains(ROI_file_list, ROIs_not_used);
    ROIfiles_lh = ROI_file_list(file_mask & contains(ROI_file_list, '_lh_')); % all lh ROIs for subj
    ROIfiles_rh = ROI_file_list(file_mask & contains(ROI_file_list, '_rh_')); % all rh ROIs for subj
    
    if isempty(ROIfiles_rh) && isempty(ROIfiles_lh)
        disp(['Did not find any ROI files for subj ' subjCode ' ... skipping ...']);
        continue
    end

    %% lh annot file creation
    % Initialize variables necessary for creating lh annot file
    annot_labels_lh = zeros(length(verts_ref),1); % start with labels as all 0s
    annot_ctable_lh.numEntries = length(all_ROIs_use) + 1; % 12 ROIs plus one for 'unknown' label 0
    annot_ctable_lh.orig_tab = ctable_ref.orig_tab; % keep same colortable reference
    annot_ctable_lh.struct_names = {ctable_ref.struct_names{1}}; % 'unknown' label is always 1st
    annot_ctable_lh.table = ctable_ref.table(1,:); % copy canonical color RGB and label for 'unknown' label
    hasROIs = cell(length(ROIfiles_lh),1);
    n_rois_lh = length(ROIfiles_lh);

    % Loop through lh ROIs
    for rr = 1:n_rois_lh

        % Load ROI
        ROI_file = ROIfiles_lh{rr};
        ROI_data_lh = MRIread([projectDir '/data/ROIs/' ROI_file]); % read ROI data
        
        % Add struct_name (label name)
        whichROI = cellfun(@(pattern) contains(ROI_file,pattern), all_ROIs_use); % determine which ROI name is in this file
        ROI_name = all_ROIs_use{whichROI};
        hasROIs{rr} = ROI_name;
        annot_ctable_lh.struct_names{rr+1} = ROI_name; % use this name for struct name (rr+1 because we added unknown struct name already

        % Add ctable row and update label index to match
        annot_ctable_lh.table(rr+1,:) = ctable_ref.table(rr+1,:); % just keep same canonical label number and color for this ROI
        label_curr = annot_ctable_lh.table(rr+1,5); % Extract the specific label number
        annot_labels_lh(logical(ROI_data_lh.vol)) = label_curr; % replace 1s with label number in annotation labels variable

    end

    % Add dummy ROIs for missing ROIs so that Conn will be able to use all ROIs
    missing_ROIs = all_ROIs_use(~ismember(all_ROIs_use, hasROIs'));
    for mr = 1:length(missing_ROIs)
        missing_ROIs_allsubj{ss}{end+1} = [missing_ROIs{mr} ' (L)'];
        annot_ctable_lh.struct_names{n_rois_lh+mr+1} = missing_ROIs{mr}; % +1 for 'unknown' name
        annot_ctable_lh.table(n_rois_lh+mr+1,:) = ctable_ref.table(n_rois_lh+mr+1,:); % just keep same canonical label number and color for this ROI
        label_curr = annot_ctable_lh.table(n_rois_lh+mr+1,5); % Extract the specific label number
        dummy_vol = false(1,fs_number); 
        dummy_vol(mr) = true; % just one voxel of ROI for dummy value
        annot_labels_lh(dummy_vol) = label_curr; % replace 1s with label number in annotation labels variable
    end

    % Save annot file
    write_annotation([projectDir '/data/ROIs/lh.' subjCode '_all_ROIs_nomissing.annot'], annot_verts, annot_labels_lh, annot_ctable_lh);

    %% rh annot file creation
    % Initialize variables necessary for creating rh annot file
    annot_labels_rh = zeros(length(verts_ref),1);
    annot_ctable_rh.numEntries = length(all_ROIs_use) + 1; % plus one for 'unknown' label 0
    annot_ctable_rh.orig_tab = ctable_ref.orig_tab;
    annot_ctable_rh.struct_names = {ctable_ref.struct_names{1}}; % 'unknown' label is always 1st
    annot_ctable_rh.table = ctable_ref.table(1,:); % copy canonical color RGB and label fpr 'unknown' label
    hasROIs = cell(length(ROIfiles_rh),1);
    n_rois_rh = length(ROIfiles_rh);

    % Loop through rh ROIs
    for rr = 1:length(ROIfiles_rh)

        % Load ROI
        ROI_file = ROIfiles_rh{rr};
        ROI_data_rh = MRIread([projectDir '/data/ROIs/' ROI_file]);
        
        % Add struct_name (label name)
        whichROI = cellfun(@(pattern) contains(ROI_file,pattern), all_ROIs_use);
        ROI_name = all_ROIs_use{whichROI};
        hasROIs{rr} = ROI_name;
        annot_ctable_rh.struct_names{rr+1} = ROI_name;

        % Add ctable row and update label index to match
        annot_ctable_rh.table(rr+1,:) = ctable_ref.table(rr+1,:);
        label_curr = annot_ctable_rh.table(rr+1,5);
        annot_labels_rh(logical(ROI_data_rh.vol)) = label_curr;

    end

    % Add dummy ROIs for missing ROIs so that Conn will be able to use all ROIs
    missing_ROIs = all_ROIs_use(~ismember(all_ROIs_use, hasROIs'));
    for mr = 1:length(missing_ROIs)
        missing_ROIs_allsubj{ss}{end+1} = [missing_ROIs{mr} ' (R)'];
        annot_ctable_rh.struct_names{n_rois_rh+mr+1} = missing_ROIs{mr}; % +1 for 'unknown' name
        annot_ctable_rh.table(n_rois_rh+mr+1,:) = ctable_ref.table(n_rois_rh+mr+1,:); % just keep same canonical label number and color for this ROI
        label_curr = annot_ctable_rh.table(n_rois_rh+mr+1,5); % Extract the specific label number
        dummy_vol = false(1,fs_number); 
        dummy_vol(1) = true; % just one voxel of ROI for dummy value
        annot_labels_rh(dummy_vol) = label_curr; % replace 1s with label number in annotation labels variable
    end

    % Save annot file
    write_annotation([projectDir '/data/ROIs/rh.' subjCode '_all_ROIs_nomissing.annot'], annot_verts, annot_labels_rh, annot_ctable_rh);

    %% Finished
    disp(['Finished creating ROI annot files for subj ' subjCode]);
end

