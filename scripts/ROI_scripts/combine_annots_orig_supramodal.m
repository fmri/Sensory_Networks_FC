%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to combine the labels for auditory biased,
% visual biased, supramodal, and posterior auditory/visual ROIs in the same
% annotation file with careful overlap control to maintain >100 vertex ROIs
%
% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Initialize parameters
projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'}; % first 20 are localizer analysis subjs, AI is in rs only
N = length(subjCodes);
av_ROIs = {'sPCS', 'iPCS', 'midIFS', 'pVis', 'tgPCS', 'cIFSG', 'CO', 'FO', 'pAud'};
sm_ROIs = {'sm_sPCS', 'sm_iPCS', 'sm_midFSG', 'sm_aINS','sm_preSMA', 'sm_dACC'};
N_ROIs = length(sm_ROIs);
hemis = {'lh', 'rh'};
fs_number = 163842; % vertices in a hemisphere

%% Read in reference annot file so we can keep the same structure for our annot files
annotpath = [projectDir 'data/recons/fsaverage/label/lh.aparc.annot'];
[verts_ref, labels_ref, ctable_ref] = read_annotation(annotpath); % Read in lh.aparc.annot file for reference annotation file structure

%% Loop through subjects and combine labels into single annotation file
for ss = 1:N
    subjCode = subjCodes{ss};
    for hh = 1:2
        hemi = hemis{hh};
        % Load subj annotation files for sm and av ROIs
        [~, labels_orig, ctable_orig] = read_annotation([ROI_dir hemi '.' subjCode '_ROIs.annot']); % We will add to this annot variable structure
        [~, labels_sm, ctable_sm] = read_annotation([ROI_dir hemi '.' subjCode '_smROIs2.annot']); 
        ctable_orig.numEntries = ctable_orig.numEntries + length(sm_ROIs);

        for rr = 1:N_ROIs

            sm_ROI = sm_ROIs{rr};
            sm_ROI_mask = labels_sm == ctable_sm.table(strcmp(ctable_sm.struct_names, sm_ROI),5); % Get vertices for this ROI 

            % Add struct_name (label name)
            curr_ROI_row = height(ctable_orig.table);
            ctable_orig.struct_names{curr_ROI_row+1} = sm_ROI; % use this name for struct name (rr+1 because we added unknown struct name already)

            % Add ctable row and update label index to match
            ctable_orig.table(curr_ROI_row+1,:) = ctable_ref.table(curr_ROI_row+1,:); % just keep same canonical label number and color for this ROI
            label_curr = ctable_orig.table(curr_ROI_row+1,5); % Extract the specific label number

            % Make sure there is no overlap
            prelabeled_verts = labels_orig(sm_ROI_mask);
            already_labeled = prelabeled_verts ~= 0; % find vertices in this ROI that already have a label
            already_labeled_vals = prelabeled_verts(already_labeled); % get the label value for verts that have label already
            if ~isempty(already_labeled_vals)
                overlapped_ROI_labels = unique(already_labeled_vals);
                N_overlapped = length(overlapped_ROI_labels);

                % create overlap safe mask (mask that ensures all previous ROIs do not get overlapped by current ROI)
                overlap_safe_mask = sm_ROI_mask;
                curr_nverts = sum(sm_ROI_mask);
                reduced_curr_lab_nverts = curr_nverts;
                for ii = 1:N_overlapped
                    overlap_safe_mask = overlap_safe_mask & labels_orig~=overlapped_ROI_labels(ii);
                    reduced_curr_lab_nverts = reduced_curr_lab_nverts - sum(already_labeled_vals==overlapped_ROI_labels(ii));
                end

                % Loop through previous ROIs that may get overlapped by current ROI and decide which way to overlap
                for ii = 1:N_overlapped
                    overlapped_ROI_name = ctable_orig.struct_names(ctable_orig.table(:,5)==overlapped_ROI_labels(ii));
                    overlapping_label_nverts = sum(labels_orig==overlapped_ROI_labels(ii)); 
                    numoverlap = sum(already_labeled_vals==overlapped_ROI_labels(ii));
                    reduced_prev_lab_nverts = overlapping_label_nverts - numoverlap;
                    if ismember(sm_ROI, {'sm_aINS', 'sm_preSMA', 'sm_dACC'}) && ismember(overlapped_ROI_name, {'aINS', 'preSMA', 'dACC'})
                        labels_orig(sm_ROI_mask) = label_curr; % we always want the full sm version of these
                        if ii ~= N_overlapped
                            overlap_safe_mask = sm_ROI_mask; % redo the overlap safe mask etc without the aINs/preSMA/dACC
                            reduced_curr_lab_nverts = curr_nverts;
                            overlap_safe_mask = overlap_safe_mask & ismember(labels_orig, overlapped_ROI_labels(ii+1:end));
                            reduced_curr_lab_nverts = reduced_curr_lab_nverts - sum(ismember(already_labeled_vals, overlapped_ROI_labels(ii+1:end)));
                        end
                    elseif reduced_prev_lab_nverts < 100 && reduced_curr_lab_nverts < 100
                        if strcmp(subjCode, 'TP') && strcmp(hemi, 'lh') && strcmp(sm_ROI, 'sm_midFSG') && strcmp(overlapped_ROI_name{:}, 'midIFS')
                            labels_orig(overlap_safe_mask) = label_curr;
                        else
                            error(['Subj ' subjCode ' ' hemi ' ' sm_ROI ' overlaps with ' overlapped_ROI_name{:} ' which causes one ROI to be less than 100 vertices. This should not happen, check label files for both ROIs.'])
                        end
                    elseif reduced_prev_lab_nverts < reduced_curr_lab_nverts % don't relabel overlap if it would leave the other ROI with fewer vertices than this one
                        labels_orig(overlap_safe_mask) = label_curr;
                    elseif reduced_curr_lab_nverts <= reduced_prev_lab_nverts % do relabel overlap if this ROI is smaller
                        labels_orig(sm_ROI_mask) = label_curr;
                    else
                        error('Encountered unexpected condition')
                    end
                end
            else
                labels_orig(sm_ROI_mask) = label_curr; % replace all indices with current ROI label number
            end
        end

        % Actually write annotation variable to file
        write_annotation([projectDir '/data/ROIs/' hemi '.' subjCode '_avsm_ROIs.annot'], verts_ref, labels_orig, ctable_orig);

    end
end







