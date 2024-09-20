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

ROIs = {'aINS', 'preSMA', 'ppreCun', 'dACC', ... % multisensory
    'sPCS', 'iPCS', 'midIFS', 'pVis_mod' ... % visual
    'tgPCS', 'cIFSG', 'pAud', 'CO', 'FO', 'cmSFG'}; % auditory
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
        subj_labels = all_label_files(contains(all_label_files, [subjCode '_']) & ~contains(all_label_files, '_indiv') ...
            & contains(all_label_files, ['_' hemi]) & contains(all_label_files, ROIs) & contains(all_label_files, '.label'));
        if strcmp(subjCode, 'NS') % due to "NS" being in the label aINS, we have to take out the wrong aINS for this subj
            bad_aINS = contains(subj_labels, 'aINS') & ~contains(subj_labels, 'NS_aINS');
            subj_labels = subj_labels(~bad_aINS);
        end

        % Check for replacement ROIs and make sure they are used
        replacement_ROIs = contains(subj_labels, 'replacement.label');
        if any(replacement_ROIs)
            rep_ROI_inds = find(replacement_ROIs);
            bad_nonrepl = false(length(subj_labels),1);
            for ii = 1:length(rep_ROI_inds)
                repl_ROI_name = strsplit(subj_labels{rep_ROI_inds(ii)}, '_');
                repl_ROI_name = repl_ROI_name{2};
                bad_nonrepl(contains(subj_labels, '.nii') | contains(subj_labels, [repl_ROI_name '_' hemi '.label'])) = true;
            end
            subj_labels = subj_labels(~bad_nonrepl);
        end

        % Check for missing ROIs
        any_missing = false;
        if length(subj_labels) ~= N_ROIs
            any_missing = true;
            subj_label_namesplit = cellfun(@(x) strsplit(x,'_'), subj_labels,'UniformOutput',false);
            subj_label_names = cellfun(@(x) x{2}, subj_label_namesplit, 'UniformOutput',false);
            missing_ROIs = ROIs(~contains(ROIs, subj_label_names));
            disp(['Subj ' subjCode ' missing ROI ' missing_ROIs{:} ' ' hemi ])
        end

        % Initialize variables necessary for creating lh annot file
        annot_labels = zeros(length(verts_ref),1); % start with labels as all 0s
        annot_ctable.numEntries = length(ROIs) + 1; % ROIs plus one for 'unknown' label 0
        annot_ctable.orig_tab = ctable_ref.orig_tab; % keep same colortable reference
        annot_ctable.struct_names = {ctable_ref.struct_names{1}}; % 'unknown' label is always 1st
        annot_ctable.table = ctable_ref.table(1,:); % copy canonical color RGB and label for 'unknown' label

        % Loop through ROIs
        N_ROIs_subj = length(subj_labels);
        for rr = 1:N_ROIs_subj
            % Load ROI
            ROI_file = subj_labels{rr};
            label_data = readtable([projectDir '/data/ROIs/' ROI_file], 'FileType','text'); % read ROI data
            label_vertex_inds = label_data{:,1} + 1; % inds off by one
            if length(label_vertex_inds) < 100
                keyboard;
            end

            % Add struct_name (label name)
            whichROI = cellfun(@(x) contains(ROI_file,x), ROIs); % determine which ROI name is in this file
            ROI_name = ROIs{whichROI};
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
                    if reduced_prev_lab_nverts < 100 && reduced_curr_lab_nverts < 100
                        error(['Subj ' subjCode ' ' hemi ' ' ROI_name ' overlaps with a label that causes both ROIs to be less than 100 vertices. This should not happen, check label files for both ROIs.'])
                    elseif reduced_prev_lab_nverts < 100 % don't relabel overlap if it will leave the other ROI with < 100 verts
                        annot_labels(overlap_safe_mask) = label_curr;
                    elseif reduced_curr_lab_nverts < 100 || ~ismember(ROI_name, {'aINS', 'preSMA', 'ppreCun', 'dACC'}) % do relabel overlap if not doing so willleave the current ROI with < 100 verts ... or if both will survive above 100 verts and current label is not multisensory, replace overlap with current ROI
                        annot_labels(label_vertex_inds) = label_curr; % replace all indices with current ROI label number
                    elseif ismember(ROI_name, {'aINS', 'preSMA', 'ppreCun', 'dACC'}) % if current ROI is multisensory, don't relabel overlap
                        annot_labels(overlap_safe_mask) = label_curr;
                    else
                        error('Encountered unexpected condition')
                    end
                end
            else
                annot_labels(label_vertex_inds) = label_curr; % replace all indices with current ROI label number
            end
        end

        % If there are any missing ROIs, add a dummy ROI to the annotation
        if any_missing
            for mr = 1:length(missing_ROIs)
                annot_ctable.struct_names{N_ROIs_subj+mr+1} = missing_ROIs{mr}; % +1 for 'unknown' name
                annot_ctable.table(N_ROIs_subj+mr+1,:) = ctable_ref.table(N_ROIs_subj+mr+1,:); % just keep same canonical label number and color for this ROI
                label_curr = annot_ctable.table(N_ROIs_subj+mr+1,5); % Extract the specific label number
                dummy_vol = false(1,fs_number);
                dummy_vol(mr) = true; % just one voxel of ROI for dummy value
                annot_labels(dummy_vol) = label_curr; % replace 1s with label number in annotation labels variable
            end
        end

        % if ~contains('VOT', annot_ctable.struct_names)
        %     keyboard;
        % end

        % VOT_labelnum = annot_ctable.table(contains(annot_ctable.struct_names, 'VOT'),5);
        % num_verts_VOT = sum(annot_labels == VOT_labelnum);
        % disp(num_verts_VOT)
        % if num_verts_VOT < 100
        %     keyboard;
        % end

        % Actually write annotation file
        write_annotation([projectDir '/data/ROIs/' hemi '.' subjCode '_ROIs_pVisMod.annot'], annot_verts, annot_labels, annot_ctable);

    end

    disp(['Finished subj ' subjCode]);

end









