%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to combine the labels for auditory biased,
% visual biased, supramodal, and posterior auditory/visual ROIs in the same
% annotation file with careful overlap control to maintain >100 vertex ROIs
%
% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/'));
ccc;

%% Initialize parameters
projectDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/';
ROI_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs/';
subj_data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'}; 
N = length(subjCodes);
av_ROIs = {'sPCS', 'iPCS', 'midIFS', 'pVis', 'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG', 'pAud'};
sm_ROIs = {'sm_sPCS', 'sm_iPCS', 'sm_midFSG', 'sm_aINS','sm_preSMA', 'sm_dACC'};
N_ROIs = length(sm_ROIs);
hemis = {'lh', 'rh'};
t_thresh = 2;
fs_number = 163842; % vertices in a hemisphere
badlist = {};
test1 = [];

%% Load supramodal search spaces
probROI_data = cell(N_ROIs, 2);
for rr = 1:N_ROIs
    for hh = 1:2
        path = [ROI_dir hemis{hh} '_intersect_' sm_ROIs{rr}(4:end) '_supra_vA-vP+aA-aP_probabilistic_thresh5_dilate5.label'];
        probROI_data{rr,hh} = readtable(path, 'FileType','text');
    end
end


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
        [~, labels_sm, ctable_sm] = read_annotation([ROI_dir hemi '.' subjCode '_smROIs.annot']); 
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

            % Deal with overlap
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
                    
                    disp(['Subj ' subjCode ' ' hemi ' ' sm_ROI ' (' num2str(curr_nverts) ') overlaps ' overlapped_ROI_name{:} ' (' num2str(overlapping_label_nverts) ') by ' num2str(numoverlap) ' verts'])

                    if ismember(overlapped_ROI_name, {'aINS', 'preSMA', 'dACC'}) % these ROIs are not actually used, so always overlap them
                        if ii ~= N_overlapped
                            overlap_safe_mask = sm_ROI_mask; % redo the overlap safe mask etc without the aINs/preSMA/dACC
                            reduced_curr_lab_nverts = curr_nverts;
                            overlap_safe_mask = overlap_safe_mask & ~ismember(labels_orig, overlapped_ROI_labels(ii+1:end));
                            reduced_curr_lab_nverts = reduced_curr_lab_nverts - sum(ismember(already_labeled_vals, overlapped_ROI_labels(ii+1:end)));
                            final_sm_ROI_mask = overlap_safe_mask;
                            labels_orig(overlap_safe_mask) = label_curr; 
                        else
                            final_sm_ROI_mask = sm_ROI_mask;
                            labels_orig(sm_ROI_mask) = label_curr; 
                        end
                    elseif reduced_prev_lab_nverts < 100 && reduced_curr_lab_nverts < 100
                        if strcmp(subjCode, 'TP') && strcmp(hemi, 'lh') && strcmp(sm_ROI, 'sm_midFSG') && strcmp(overlapped_ROI_name{:}, 'midIFS')
                            final_sm_ROI_mask = overlap_safe_mask;
                            labels_orig(overlap_safe_mask) = label_curr;
                            continue;
                        end
                        error(['Subj ' subjCode ' ' hemi ' ' sm_ROI ' overlaps ' overlapped_ROI_name{:} ' which causes both ROIs to be less than 100 vertices. This should not happen, check label files for both ROIs.'])
                    elseif reduced_prev_lab_nverts < 100 % if the sensory-biased ROI will be less than 100 vertices if overlapped, give all vertices to sensory-biased ROI
                        final_sm_ROI_mask = overlap_safe_mask;
                        labels_orig(overlap_safe_mask) = label_curr;
                    elseif reduced_curr_lab_nverts <= 100 % If supramodal ROI will be less than 100 vertices if overlapped, still exclude sensory-biased vertices and find the next strongest supramodal vertices to add to supramodal ROI so that it is above 100
                        verts_needed = 100 - reduced_curr_lab_nverts;
                        probROI = probROI_data{rr,hh};
                        prev_label_mask = labels_orig==overlapped_ROI_labels(ii);
                        search_space = ismember( 1:fs_number, probROI.Var1+1);
                        search_space_lim = search_space' & ~prev_label_mask & ~sm_ROI_mask;
                        search_space_lim_inds = find(search_space_lim);
                        path_vA_vP = [subj_data_dir subjCode '/localizer/localizer_contrasts_' hemi '/vA-vP/t.nii.gz'];
                        path_aA_aP = [subj_data_dir subjCode '/localizer/localizer_contrasts_' hemi '/aA-aP/t.nii.gz'];
                        path_V_A = [subj_data_dir subjCode '/localizer/localizer_contrasts_' hemi '/V-A/t.nii.gz'];
                        vA_vP = MRIread(path_vA_vP);
                        aA_aP = MRIread(path_aA_aP);
                        V_A = MRIread(path_V_A);
                        vA_vP_inROI = vA_vP.vol(search_space_lim); 
                        aA_aP_inROI = aA_aP.vol(search_space_lim); 
                        V_A_inROI = V_A.vol(search_space_lim); 
                        candidate_verts = abs(V_A_inROI) < 2 & vA_vP_inROI > t_thresh & aA_aP_inROI > t_thresh;
                        if sum(candidate_verts) < verts_needed
                            badlist{end+1} =  [subjCode '_' sm_ROI '_' hemi ];
                            t_thresh2 = t_thresh;
                            while sum(candidate_verts) < verts_needed
                                t_thresh2 = t_thresh2 - 0.001;
                                candidate_verts = abs(V_A_inROI) < 2 & vA_vP_inROI > t_thresh2 & aA_aP_inROI > t_thresh2;
                            end
                            disp(['t thresh needed to get 100 verts for supramodal ROI = ' num2str(t_thresh2)]);
                        end
                        verts_add = ismember(1:fs_number, search_space_lim_inds(candidate_verts));
                        final_sm_ROI_mask = verts_add' | overlap_safe_mask;
                        labels_orig(final_sm_ROI_mask) = label_curr;
                    else % if neither will dip below 100 vertices, give the verts to sensory-biased ROI bc one modality is more responsive
                        final_sm_ROI_mask = overlap_safe_mask;
                        labels_orig(overlap_safe_mask) = label_curr;
                    end
                end
            else
                final_sm_ROI_mask = sm_ROI_mask;
                labels_orig(sm_ROI_mask) = label_curr; % replace all indices with current ROI label number
            end
        end

        % Actually write annotation variable to file
        write_annotation([projectDir '/data/ROIs/' hemi '.' subjCode '_avsm_ROIs.annot'], verts_ref, labels_orig, ctable_orig);

    end
end







