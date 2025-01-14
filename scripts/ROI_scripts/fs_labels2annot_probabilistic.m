%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create freesurfer annotation files
%%% from the probabilistic label files and eliminate overlap between ROIs
%%% Tom Possidente - December 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;



%% Set key variables
% Read in reference annot file so we can keep the same structure for our annot files
annotpath = '/projectnb/somerslab/tom/projects/spacetime_network/data/recons/fsaverage/label/lh.aparc.annot';
[verts_ref, labels_ref, ctable_ref] = read_annotation(annotpath); % Read in lh.aparc.annot file for reference annotation file structure

% Create vertices for all annot files to use (same number of vertices for all bc fsaverage space)
annot_verts = verts_ref;

ROIs = {'sPCS', 'iPCS', 'midIFS', 'pVis' ... % visual
    'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG', 'pAud',... % auditory
    'aINS', 'dACC', 'preSMA'}; % multisensory 
N_ROIs = length(ROIs);

fs_number = 163842;
ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
hemis = {'lh', 'rh'};

%% Load in probabilistic ROIs
prob_ROI_data = cell(N_ROIs,2);
for rr = 1:N_ROIs
    ROI = ROIs{rr};
    for hh = 1:2
        thresh = 5; % N=5 threshold on probabilistic ROIs
        prob_ROI_fpath = [ROI_dir 'probabilistic_' ROI '_' hemis{hh} '_thresh' num2str(thresh) '_dilated5.label'];
        while ~isfile(prob_ROI_fpath)
            thresh = thresh - 1;
            prob_ROI_fpath = [ROI_dir 'probabilistic_' ROI '_' hemis{hh} '_thresh' num2str(thresh) '_dilated5.label'];
            if thresh < 3
                error(['probabilistic ROI for ' ROI ' not found'])
            end
        end
        prob_ROI_data{rr,hh} = readtable(prob_ROI_fpath, 'FileType', 'text');
    end
end

%% create annot file for each hemisphere, while eliminating overlap
for hh = 1:length(hemis)

    % Initialize variables necessary for creating lh annot file
    annot_labels = zeros(length(verts_ref),1); % start with labels as all 0s
    annot_ctable.numEntries = length(ROIs) + 1; % ROIs plus one for 'unknown' label 0
    annot_ctable.orig_tab = ctable_ref.orig_tab; % keep same colortable reference
    annot_ctable.struct_names = {ctable_ref.struct_names{1}}; % 'unknown' label is always 1st
    annot_ctable.table = ctable_ref.table(1,:); % copy canonical color RGB and label for 'unknown' label

    for rr = 1:N_ROIs
        % Load ROI
        ROI_data = prob_ROI_data{rr,hh};
        label_vertex_inds = ROI_data.Var1+1;

        % Add struct_name (label name)
        ROI_name = ROIs{rr};
        annot_ctable.struct_names{rr+1} = ROI_name; % use this name for struct name (rr+1 because we added unknown struct name already)

        % Add ctable row and update label index to match
        annot_ctable.table(rr+1,:) = ctable_ref.table(rr+1,:); % just keep same canonical label number and color for this ROI
        label_curr = annot_ctable.table(rr+1,5); % Extract the specific label number

        % Deal with possible overlap
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
                if reduced_prev_lab_nverts < 100 || reduced_curr_lab_nverts < 100
                    error(['Subj ' subjCode ' ' hemi ' ' ROI_name ' overlaps with a label that causes one ROI to be less than 100 vertices. This should not happen, check label files for both ROIs.'])
                elseif ismember(ROI_name, {'aINS', 'preSMA', 'ppreCun', 'dACC'}) % if current ROI is multisensory, don't relabel overlap
                    annot_labels(overlap_safe_mask) = label_curr;
                elseif reduced_prev_lab_nverts < reduced_curr_lab_nverts % don't relabel overlap if it would leave the other ROI with fewer vertices than this one
                    annot_labels(overlap_safe_mask) = label_curr;
                elseif reduced_curr_lab_nverts < reduced_prev_lab_nverts % do relabel overlap if this ROI is smaller
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
    write_annotation([ROI_dir hemis{hh} '.probabilistic_ROIs.annot'], annot_verts, annot_labels, annot_ctable);
end











