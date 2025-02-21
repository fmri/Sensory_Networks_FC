%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create freesurfer annotation files
%%% from the freesurfer label files for each subj
%%% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Get Subj Codes

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
experiment_name = 'spacetime';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
reject_subjs = {'RR', 'AH', 'SL'};
subjCodes = subjCodes(~ismember(subjCodes, reject_subjs));
N = length(subjCodes);

%% Set key variables
% Read in reference annot file so we can keep the same structure for our annot files
annotpath = [projectDir 'data/recons/fsaverage/label/lh.aparc.annot'];
[verts_ref, labels_ref, ctable_ref] = read_annotation(annotpath); % Read in lh.aparc.annot file for reference annotation file structure

% Create vertices for all annot files to use (same number of vertices for all bc fsaverage space)
annot_verts = verts_ref;

ROIs = {'sm_sPCS', 'sm_iPCS', 'sm_midFSG', 'sm_aINS', 'sm_preSMA', 'sm_dACC', 'pVis', 'pAud'};
N_ROIs = length(ROIs);

fs_number = 163842;

%% Loop through subjs and create annot file
ROI_dir = [projectDir '/data/ROIs/'];

all_label_files = {dir(ROI_dir).name};
hemis = {'lh', 'rh'};

for ss = 1:N

    subjCode = subjCodes{ss};
    for hh = 1:length(hemis)

        % Get labels for this subj/hemisphere
        hemi = hemis{hh};
        subj_labels = all_label_files(contains(all_label_files, [subjCode '_' hemi '_sm_']) & contains(all_label_files, ROIs) & contains(all_label_files, '2.label'));
        subj_posterior_labels = {[subjCode '_pVis_' hemi '.label'],[subjCode '_pAud_' hemi '.label']};
        subj_labels = cat(2, subj_labels, subj_posterior_labels);

        % Check for missing ROIs
        assert(length(subj_labels) == N_ROIs, 'Number of subj label files does not match number of ROIs')

        % Initialize variables necessary for creating annot file
        annot_labels = zeros(length(verts_ref),1); % start with labels as all 0s
        annot_ctable.numEntries = length(subj_labels) + 1; % ROIs plus one for 'unknown' label 0
        annot_ctable.orig_tab = ctable_ref.orig_tab; % keep same colortable reference
        annot_ctable.struct_names = {ctable_ref.struct_names{1}}; % 'unknown' label is always 1st
        annot_ctable.table = ctable_ref.table(1,:); % copy canonical color RGB and label for 'unknown' label

        % Loop through ROIs
        for rr = 1:length(subj_labels)

            % Load ROI
            ROI_file = subj_labels{rr};
            label_data = readtable([projectDir '/data/ROIs/' ROI_file], 'FileType','text'); % read ROI data
            label_vertex_inds = label_data{:,1} + 1; % inds off by one
            assert(length(label_vertex_inds) >= 100, 'ROI has fewer than 100 vertices')

            % Add struct_name (label name)
            whichROI = cellfun(@(x) contains(ROI_file,x), ROIs); % determine which ROI name is in this file
            ROI_name = ROIs{whichROI};
            annot_ctable.struct_names{rr+1} = ROI_name; % use this name for struct name (rr+1 because we added unknown struct name already)

            % Add ctable row and update label index to match
            annot_ctable.table(rr+1,:) = ctable_ref.table(rr+1,:); % just keep same canonical label number and color for this ROI
            label_curr = annot_ctable.table(rr+1,5); % Extract the specific label number

            % Make sure there is no overlap
            prelabeled_verts = annot_labels(label_vertex_inds);
            already_labeled = prelabeled_verts ~= 0; % find vertices in this ROI that already have a label
            already_labeled_vals = prelabeled_verts(already_labeled); % get the label value for verts that have label already
            if ~isempty(already_labeled_vals)
                overlapped_ROI_labels = unique(already_labeled_vals);

                % create overlap safe mask (mask that ensures all previous ROIs do not get overlapped by current ROI)
                inds = 1:fs_number;
                label_vertex_mask = ismember(inds, label_vertex_inds); % convert from indices to binary mask
                overlap_safe_mask = label_vertex_mask';
                curr_nverts = length(label_vertex_inds);
                reduced_curr_lab_nverts = curr_nverts;
                for ii = 1:length(overlapped_ROI_labels)
                    overlap_safe_mask = overlap_safe_mask & annot_labels~=overlapped_ROI_labels(ii);
                    reduced_curr_lab_nverts = reduced_curr_lab_nverts - sum(already_labeled_vals==overlapped_ROI_labels(ii));
                end

                % Loop through previous ROIs that may get overlapped by current ROI and decide which way to overlap
                for ii = 1:length(overlapped_ROI_labels)
                    overlapping_label_nverts = sum(annot_labels==overlapped_ROI_labels(ii));
                    numoverlap = sum(already_labeled_vals==overlapped_ROI_labels(ii));
                    reduced_prev_lab_nverts = overlapping_label_nverts - numoverlap;
                    if reduced_prev_lab_nverts < 100 && reduced_curr_lab_nverts < 100
                        error(['Subj ' subjCode ' ' hemi ' ' ROI_name ' overlaps with a label that causes one ROI to be less than 100 vertices. This should not happen, check label files for both ROIs.'])
                    elseif reduced_prev_lab_nverts < reduced_curr_lab_nverts % don't relabel overlap if it would leave the other ROI with fewer vertices than this one
                        annot_labels(overlap_safe_mask) = label_curr;
                    elseif reduced_curr_lab_nverts <= reduced_prev_lab_nverts % do relabel overlap if this ROI is smaller
                        annot_labels(label_vertex_inds) = label_curr;
                    else
                        error('Encountered unexpected condition')
                    end
                end
            else
                annot_labels(label_vertex_inds) = label_curr; % replace all indices with current ROI label number
            end
        end

        % Actually write annotation file
        write_annotation([projectDir '/data/ROIs/' hemi '.' subjCode '_smROIs2.annot'], annot_verts, annot_labels, annot_ctable);

    end

    disp(['Finished subj ' subjCode]);

end









