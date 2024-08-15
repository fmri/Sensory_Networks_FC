%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create freesurfer annotation files
%%% from the freesurfer label files for each subj
%%% Tom Possidente - August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Get Subj Codes

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
experiment_name = 'spacetime';


subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'SL', 'AH', 'RR'})); % rejected subjs
N = length(subjCodes);

%% Set key variables
% Read in reference annot file so we can keep the same structure for our annot files
annotpath = [projectDir 'data/recons/fsaverage/label/lh.aparc.annot'];
[verts_ref, labels_ref, ctable_ref] = read_annotation(annotpath); % Read in lh.aparc.annot file for reference annotation file structure

% Create vertices for all annot files to use (same number of vertices for all bc fsaverage space)
annot_verts = verts_ref;

ROIs = ["aINS", "preSMA", "ppreCun", "dACC", ... % multisensory
        "sPCS", "iPCS", "midIFS", "aIPS", "pIPS", "DO", "LOT", "VOT",... % visual
        "tgPCS", "cIFSG", "pAud", "CO", "FO", "cmSFG"]'; % auditory
N_ROIs = length(ROIs);

%% Loop through subjs and create annot file
ROI_dir = [projectDir '/data/ROIs/'];

all_label_files = {dir(ROI_dir).name};
hemis = {'lh', 'rh'};

for ss = 1:N

    subjCode = subjCodes{ss};
    for hh = 1:length(hemis)

        % Get labels for this subj/hemisphere
        hemi = hemis{hh};
        subj_labels = all_label_files(contains(all_label_files, [subjCode '_']) & ~contains(all_label_files, '_indiv') ...
            & contains(all_label_files, ['_' hemi]));
        if strcmp(subjCode, 'NS') % due to "NS" being in the label aINS, we have to take out the wrong aINS for this subj
            bad_aINS = contains(subj_labels, 'aINS') & ~contains(subj_labels, 'NS_aINS');
            subj_labels = subj_labels(~bad_aINS);
        end

        assert(length(subj_labels)==N_ROIs, ['Number of ROIs for ' subjCode ' ' hemi ' is ' num2str(length(subj_labels)) ...
                                            ' but should be ' num2str(N_ROIs)]);

        % Initialize variables necessary for creating lh annot file
        annot_labels = zeros(length(verts_ref),1); % start with labels as all 0s
        annot_ctable.numEntries = length(ROIs) + 1; % ROIs plus one for 'unknown' label 0
        annot_ctable.orig_tab = ctable_ref.orig_tab; % keep same colortable reference
        annot_ctable.struct_names = {ctable_ref.struct_names{1}}; % 'unknown' label is always 1st
        annot_ctable.table = ctable_ref.table(1,:); % copy canonical color RGB and label for 'unknown' label

        % Loop through ROIs
        for rr = 1:N_ROIs

            % Load ROI
            ROI_file = subj_labels{rr};
            label_data = readtable([projectDir '/data/ROIs/' ROI_file], 'FileType','text'); % read ROI data
            label_vertex_inds = label_data{:,1};

            % check of there are every 0s or 163842 (if there are 0s we are
            % good, if there is a 163842 we are off by one and should
            % subtract 1 from label_vertex_inds 
            if any(ismember([0, 163842], label_vertex_inds))
                keyboard
            end

            % Add struct_name (label name)
            whichROI = cellfun(@(x) contains(ROI_file,x), ROIs); % determine which ROI name is in this file
            ROI_name = ROIs{whichROI};
            annot_ctable.struct_names{rr+1} = ROI_name; % use this name for struct name (rr+1 because we added unknown struct name already)

            % Add ctable row and update label index to match
            annot_ctable.table(rr+1,:) = ctable_ref.table(rr+1,:); % just keep same canonical label number and color for this ROI
            label_curr = annot_ctable.table(rr+1,5); % Extract the specific label number
            annot_labels(label_vertex_inds) = label_curr; % replace 0s with label number in annotation labels variable

        end
        
        % Actually write annotation file
        %write_annotation([projectDir '/data/ROIs/' hemi '.' subjCode '_ROIs.annot'], annot_verts, annot_labels, annot_ctable);
        keyboard;

    end

end









