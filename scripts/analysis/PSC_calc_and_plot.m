%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract % signal change data from
% within-subj contrasts, average them over each ROI, and finally average
% them over subjs
% Tom Possidente - September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Intialize which data we are using
data_use = 'localizer'; % Choices: 'spacetime' or 'localizer'
use_replacement_ROIs = false; 
use_probabilistic_ROIs = false; % if false, uses individualized ROIs
disp(['Calculating PSCs for ' upper(data_use) ' data']);

switch data_use
    case 'spacetime'
        task = 'spacetime';
        subjdata_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/';
        contrasts = {'sV-pV', 'tV-pV', 'sA-pA', 'tA-pA', 'pV-f', 'pA-f', 'vP-aP', 'aP-vP'};
        reject_contrast = {{}, {'RT', 'LA'}, {'LN', 'GG', 'TP'}, {}, {}, {}, {}, {}}; % These subjs had below 55% accuracy on these tasks and should be rejected from analysis in those conditions
        contrast_fpath_strs = {'bold', 'spacetime_contrasts_', '_newcondfiles'};
        if use_probabilistic_ROIs
            fullsubj_reject = {'RR', 'SL', 'AI'}; % can keep AH if doing probabilistic ROIs, AH just doesn't have enough localizer data to get good individual subj ROIs
        else
            fullsubj_reject = {'RR', 'SL', 'AI', 'AH'}; % can keep AH if doing probabilistic ROIs, AH just doesn't have enough localizer data to get good individual subj ROIs
        end
    case 'localizer'
        task = 'x1WayLocalizer';
        subjdata_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii_fs_localizer/';
        contrasts = {'vA-vP', 'aA-aP', 'f-aP', 'f-vP', 'vP-aP', 'aP-vP'};
        reject_contrast = {{}, {}, {}, {}, {}, {}, {}, {}}; % These subjs had below 55% accuracy on these tasks and should be rejected from analysis in those conditions
        contrast_fpath_strs = {'localizer', 'localizer_contrasts_', ''};
        if use_probabilistic_ROIs
            fullsubj_reject = {'RR', 'SL'};
        else
            fullsubj_reject = {'RR', 'AH', 'SL'};
        end
    otherwise 
        error('data_use variable not recognized - should be either "localizer" or "spacetime"')
end

%% Load subject info
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([task 'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, fullsubj_reject));
N = length(subjCodes);

%% Initialize key variables
ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
annot_file_suffix = '_ROIs.annot';

N_contrasts = length(contrasts);

ROIs = {'aINS', 'preSMA', 'dACC', ... % multisensory
    'sPCS', 'iPCS', 'midIFS', 'pVis' ... % visual
    'tgPCS', 'cIFSG', 'pAud', 'CO', 'FO', 'cmSFG'}; % auditory
N_ROIs = length(ROIs);
load('/projectnb/somerslab/tom/projects/spacetime_network/data/missing_ROIs.mat', 'missing_ROIs');
load('/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/replacement_ROI_list.mat', 'replacement_ROIs');

if use_probabilistic_ROIs
    annot_fpath_lh = [ROI_dir 'lh.probabilistic_ROIs.annot'];
    annot_fpath_rh = [ROI_dir 'rh.probabilistic_ROIs.annot'];
    [~, annot_labels_hemis{1}, annot_ctable_hemis{1}] = read_annotation(annot_fpath_lh);
    [~, annot_labels_hemis{2}, annot_ctable_hemis{2}] = read_annotation(annot_fpath_rh);
end

hemis = {'lh', 'rh'};

psc_results = nan(N, 2, N_contrasts, N_ROIs); % subjs x hemis x contrasts x ROIs = 20 x 2 x 4 x 13

%% Loop through subj/hemisphere/contrast/ROI and calculate % signal change in ROI

for ss = 1:N % for subj

    subjCode = subjCodes{ss};

    for hh = 1:2 % for hemisphere
        hemi = hemis{hh};

        % load subj ROIs
        if ~use_probabilistic_ROIs
            annot_fpath = [ROI_dir hemi '.' subjCode annot_file_suffix];
            [~, annot_labels, annot_ctable] = read_annotation(annot_fpath);
    
            % Make sure all ROIs exist in annot file
            assert( all(ismember(ROIs, annot_ctable.struct_names)), ['Annotation file for subj ' subjCode ' ' hemi ' does not contain all specified ROIs']);
        else
            annot_labels = annot_labels_hemis{hh};
            annot_ctable = annot_ctable_hemis{hh};
        end

        for cc = 1:N_contrasts
            contrast = contrasts{cc};

            % Check for subjs with rejected contrasts
            if ismember(subjCode, reject_contrast{cc})
                psc_results(ss, hh, cc, :) = nan;
                continue;
            end

            % load subj contrast percent signal change (psc) data
            contrast_fpath = [subjdata_dir subjCode '/' contrast_fpath_strs{1} '/' contrast_fpath_strs{2} hemi contrast_fpath_strs{3} '/' contrast '/cespct.nii.gz'];
            if ~isfile(contrast_fpath)
                contrast_fpath = [subjdata_dir subjCode '/' contrast_fpath_strs{1} '/' contrast_fpath_strs{2} hemi '_xtra_contrasts_NoyceRep/' contrast '/cespct.nii.gz'];
                if ~isfile(contrast_fpath)
                    disp(['Contrast file not found for subj ' subjCode ' contrast ' contrast ' - SKIPPING'])
                    continue
                end
            end
            psc_data = MRIread(contrast_fpath);

            for roi = 1:N_ROIs % for ROI
                
                ROI = ROIs{roi};

                if ismember([subjCode '_' ROI '_' hemi], missing_ROIs)
                    psc_results(ss, hh, cc, roi) = nan;
                    continue;
                end

                if ~use_replacement_ROIs && ismember([subjCode '_' ROI '_' hemi], replacement_ROIs)
                    psc_results(ss,hh,cc,roi) = nan;
                    continue;
                end

                ROI_index = find(ismember(annot_ctable.struct_names, ROI));
                ROI_label = annot_ctable.table(ROI_index,5);
                vertex_mask = ismember(annot_labels, ROI_label); % We can use ismember here because annotation file vertex number order is that same as those in psc_data (nii.gz)
                assert(sum(vertex_mask)>=100, ['Subj ' subjCode ' ' hemi ' ' ROI ' is below 100 vertices']);

                % Get mean psc data from ROI mask
                if use_probabilistic_ROIs && ismember(ROI, {'aINS', 'pAud'}) % Some ROIs probabilistics are expanded outside of the non-nan borders 
                     psc_results(ss,hh,cc,roi) = mean(psc_data.vol(vertex_mask), 'all', 'omitnan'); 
                else
                    psc_results(ss,hh,cc,roi) = mean(psc_data.vol(vertex_mask)); 
                end
                assert(~isnan(psc_results(ss,hh,cc,roi)), ['PSC for subj ' subjCode ' ' hemi ' ' ROI ' computed as NaN'])

            end % end ROI
        end % end contrast
    end % end hemisphere
    disp(['Finished subj: ' subjCode])
end % end subj


%% Average results over subjs
if strcmp(data_use,'localizer')
    psc_results(:,:,3:4,:) = -psc_results(:,:,3:4,:); % reverse contrast for f-vP and f-aP 
end
save('PSC_results_localizer_no_replacements.mat', 'psc_results', 'subjCodes', 'ROIs', 'contrasts');
psc_subjavg = squeeze(mean(psc_results, 1, 'omitnan'));
psc_subjhemiavg = squeeze(mean(psc_results, [1,2], 'omitnan'));
psc_subjhemiavg_tbl = array2table(psc_subjhemiavg', 'VariableNames', contrasts, 'RowNames', ROIs);
figure; heatmap(contrasts, ROIs, table2array(psc_subjhemiavg_tbl), 'Colormap', turbo)

%% Average spatial, temporal, visual, and auditory contrasts
% visual_pscs = squeeze(mean(psc_results(:,:,[1,2],:), 3, 'omitnan'));
% auditory_pscs = squeeze(mean(psc_results(:,:,[3,4],:), 3, 'omitnan'));
% spatial_pscs = squeeze(mean(psc_results(:,:,[1,3],:), 3, 'omitnan'));
% temporal_pscs = squeeze(mean(psc_results(:,:,[2,4],:), 3, 'omitnan'));
% 
% visual_pscs_subjavg = squeeze(mean(visual_pscs, 1, 'omitnan'));
% auditory_pscs_subjavg = squeeze(mean(auditory_pscs, 1, 'omitnan'));
% spatial_pscs_subjavg = squeeze(mean(spatial_pscs, 1, 'omitnan'));
% temporal_pscs_subjavg = squeeze(mean(temporal_pscs, 1, 'omitnan'));
% 
% visual_pscs_subjhemiavg = squeeze(mean(visual_pscs_subjavg, 1, 'omitnan'));
% auditory_pscs_subjhemiavg = squeeze(mean(auditory_pscs_subjavg, 1, 'omitnan'));
% spatial_pscs_subjhemiavg = squeeze(mean(spatial_pscs_subjavg, 1, 'omitnan'));
% temporal_pscs_subjhemiavg = squeeze(mean(temporal_pscs_subjavg, 1, 'omitnan'));
% 
% %% Plot PSC per ROI/hemi
% set(groot, 'DefaultAxesTickLabelInterpreter', 'none');
% 
% % subj, hemisphere, ROI PSCs
% lh_ROIs = cellfun(@(x) [x '_lh'], ROIs, 'UniformOutput', false);
% rh_ROIs = cellfun(@(x) [x '_rh'], ROIs, 'UniformOutput', false);
% hemi_ROIs = [lh_ROIs; rh_ROIs];
% hemi_ROIs = hemi_ROIs(:)'; % interleave lh and rh
% max_psc = max(psc_results,[],'all');
% min_psc = min(psc_results,[],'all');
% 
% for ii = 1:N_contrasts
% 
%     psc_onecontrast = squeeze(psc_results(:,:,ii,:));
%     psc_flat = psc_onecontrast(:);
%     xs = repelem(1:(N_ROIs*2), N);
%     xlabs = hemi_ROIs;
% 
%     figure;
%     swarmchart(xs, psc_flat, 25, 'blue', 'filled');
%     xticks(1:(N_ROIs*2));
%     xticklabels(xlabs);
%     yline(0);
%     xlabel('ROIs');
%     ylabel('% Signal Change');
%     title(contrasts{ii});
% 
%     % Plot means too
%     hold on;
%     psc_subjmean_onecontrast = squeeze(psc_subjavg(:,ii,:));
%     pscs_subjmean_flat = psc_subjmean_onecontrast(:);
%     xs_avg = 1:(N_ROIs*2);
%     swarmchart(xs_avg, pscs_subjmean_flat, 50, 'red', 'filled');
% 
%     ylim([min_psc, max_psc]);
% 
% end
% 
% 
% 
% %% Plot ROIs on spatial-temporal and visual-auditory axes
% spatial_temporal = spatial_pscs - temporal_pscs;
% visual_auditory = visual_pscs - auditory_pscs;
% ROIs_per_modality = [4,4,6]; % MD, visual, auditory
% modality_colors = {'g', 'b', 'r'};
% modality_shapes = {'^', 'o', 'd'};
% 
% figure;
% yline(0);
% xline(0);
% hold on;
% handles = {};
% ROI_count = 1;
% for mm = 1:length(ROIs_per_modality)
%     vertex_mask = ROI_count:ROI_count+ROIs_per_modality(mm)-1;
%     modality_spatial_temporal = squeeze(spatial_temporal(:,:,vertex_mask));
%     mean_modality_st = mean(modality_spatial_temporal, 1, 'omitnan');
%     modality_visual_auditory = squeeze(visual_auditory(:,:,vertex_mask));
%     mean_modality_va = mean(modality_visual_auditory, 1, 'omitnan');
%     coords = [mean_modality_va(:), mean_modality_st(:)];
%     handles{mm} = scatter(coords(:,1), coords(:,2), 100, modality_colors{mm}, modality_shapes{mm});
%     hold on;
%     ROI_count = ROI_count + ROIs_per_modality(mm);
%     disptable = table(repelem(ROIs(vertex_mask),2)', coords(:,1), coords(:,2), 'VariableNames', {'ROI', 'v-a', 's-t'});
%     disp(disptable);
% end
% xlabel('Visual-Auditory (PSC)');
% ylabel('Spatial-Temporal (PSC)');
% legend([handles{1}, handles{2}, handles{3}], 'Multiple Demand', 'Visual', 'Auditory')
% 
% 
% %% Plot ROIs on passive visual-auditory axes
% visual_auditory_passive = squeeze(psc_results(:,:,5,:)) - squeeze(psc_results(:,:,6,:));
% ROIs_per_modality = [4,4,6]; % MD, visual, auditory
% modality_colors = {'g', 'b', 'r'};
% modality_shapes = {'^', 'o', 'd'};
% 
% figure;
% yline(0);
% xline(0);
% hold on;
% handles = {};
% ROI_count = 1;
% for mm = 1:length(ROIs_per_modality)
%     vertex_mask = ROI_count:ROI_count+ROIs_per_modality(mm)-1;
%     modality_visual_auditory_passive = visual_auditory_passive(:,:,vertex_mask);
%     mean_modality_va_passive = squeeze(mean(modality_visual_auditory_passive, 1, 'omitnan'));
%     handles{mm} = scatter(mean_modality_va_passive(:), zeros(length(mean_modality_va_passive(:)),1), 100, modality_colors{mm}, modality_shapes{mm});
%     hold on;
%     ROI_count = ROI_count + ROIs_per_modality(mm);
%     disptable = table(repelem(ROIs(vertex_mask),2)', mean_modality_va_passive(:), 'VariableNames', {'ROI', 'pV-pA'});
%     disp(disptable);
% end
% xlabel('Visual-Auditory (PSC)');
% legend([handles{1}, handles{2}, handles{3}], 'Multiple Demand', 'Visual', 'Auditory')
% 
% 
% 
% 
