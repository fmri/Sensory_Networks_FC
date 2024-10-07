%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract % signal change data from
% within-subj contrasts, average them over each ROI, and finally average
% them over subjs
% Tom Possidente - September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load subject info
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.('spacetimeRuns'),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}));
N = length(subjCodes);

%% Initialize key variables
subjdata_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/';
ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
annot_file_suffix = '_ROIs.annot';

contrasts = {'sV-pV', 'tV-pV', 'sA-pA', 'tA-pA', 'pV-f', 'pA-f'};
reject_contrast = {{}, {'RT', 'LA'}, {'LN', 'GG', 'TP'}, {}, {}, {}}; % These subjs had below 55% accuracy on these tasks and should be rejected from analysis in those conditions
N_contrasts = length(contrasts);

ROIs = {'aINS', 'preSMA', 'ppreCun', 'dACC', ... % multisensory
    'sPCS', 'iPCS', 'midIFS', 'pVis' ... % visual
    'tgPCS', 'cIFSG', 'pAud', 'CO', 'FO', 'cmSFG'}; % auditory
N_ROIs = length(ROIs);
load('/projectnb/somerslab/tom/projects/spacetime_network/data/missing_ROIs.mat', 'missing_ROIs');

hemis = {'lh', 'rh'};
Nruns = 6;
missing_runs_subjs = {'LA', 'PQ'};
good_runs = {[1,4], [1,2,4]}; % corresponds to good runs for missing_runs_subjs

psc_results = nan(N, 2, N_contrasts, N_ROIs, Nruns); % subjs x hemis x contrasts x ROIs x runs = 20 x 2 x 4 x 14 x 6

%% Loop through subj/hemisphere/contrast/ROI and calculate % signal change in ROI

for ss = 1:N % for subj

    subjCode = subjCodes{ss};
    if ~ismember(subjCode, missing_runs_subjs)
        runDir = {dir([subjdata_dir subjCode '/bold/']).name};
        nruns = sum(contains(runDir, '00'));
        run_inds = 1:nruns;
    else
        missing_mask = ismember(subjCode, missing_runs_subjs);
        run_inds = good_runs{missing_mask};
    end
    for rr = 1:length(run_inds)
        run_ind = num2str(run_inds(rr));
        for hh = 1:2 % for hemisphere
            hemi = hemis{hh};

            % load subj ROIs
            annot_fpath = [ROI_dir hemi '.' subjCode annot_file_suffix];
            [annot_verts, annot_labels, annot_ctable] = read_annotation(annot_fpath);

            % Make sure all ROIs exist in annot file
            assert( all(ismember(ROIs, annot_ctable.struct_names)), ['Annotation file for subj ' subjCode ' ' hemi ' does not contain all specified ROIs']);

            for cc = 1:N_contrasts
                contrast = contrasts{cc};

                % Check for subjs with rejected contrasts
                if ismember(subjCode, reject_contrast{cc})
                    psc_results(ss, hh, cc, :, rr) = nan;
                    continue;
                end

                % load subj contrast percent signal change (psc) data
                contrast_fpath = [subjdata_dir subjCode '/bold/spacetime_contrasts_' hemi '_newcondfiles_run' run_ind '/' contrast '/cespct.nii.gz'];
                psc_data = MRIread(contrast_fpath);

                for roi = 1:N_ROIs % for ROI

                    ROI = ROIs{roi};

                    if ismember([subjCode '_' ROI '_' hemi], missing_ROIs)
                        psc_results(ss, hh, cc, roi, rr) = nan;
                        continue;
                    end

                    ROI_index = find(ismember(annot_ctable.struct_names, ROI));
                    ROI_label = annot_ctable.table(ROI_index,5);
                    ROI_mask = annot_labels == ROI_label;
                    assert(sum(ROI_mask)>=100, ['Subj ' subjCode ' ' hemi ' ' ROI ' is below 100 vertices']);

                    % Get mean psc data from ROI mask
                    psc_results(ss,hh,cc,roi,rr) = mean(psc_data.vol(ROI_mask)); %%%% MAKE SURE THESE INDS LINE UP CORRECTLY
                    assert(~isnan(psc_results(ss,hh,cc,roi,rr)), ['PSC for subj ' subjCode ' ' hemi ' ' ROI ' computed as NaN'])

                end % end ROI
            end % end contrast
        end % end hemisphere
    end % end run
end % end subj


%% Average results over subjs
%save('PSC_results.mat', 'psc_results', 'subjCodes', 'ROIs', 'contrasts');
psc_results = squeeze(mean(psc_results, 5, 'omitnan'));
psc_subjavg = squeeze(mean(psc_results, 1, 'omitnan'));
psc_subjhemiavg = squeeze(mean(psc_results, [1,2], 'omitnan'));
psc_subjhemiavg_tbl = array2table(psc_subjhemiavg', 'VariableNames', contrasts, 'RowNames', ROIs);
figure; heatmap(contrasts, ROIs, table2array(psc_subjhemiavg_tbl), 'Colormap', turbo)

%% Average spatial, temporal, visual, and auditory contrasts
visual_pscs = squeeze(mean(psc_results(:,:,[1,2],:), 3, 'omitnan'));
auditory_pscs = squeeze(mean(psc_results(:,:,[3,4],:), 3, 'omitnan'));
spatial_pscs = squeeze(mean(psc_results(:,:,[1,3],:), 3, 'omitnan'));
temporal_pscs = squeeze(mean(psc_results(:,:,[2,4],:), 3, 'omitnan'));

visual_pscs_subjavg = squeeze(mean(visual_pscs, 1, 'omitnan'));
auditory_pscs_subjavg = squeeze(mean(auditory_pscs, 1, 'omitnan'));
spatial_pscs_subjavg = squeeze(mean(spatial_pscs, 1, 'omitnan'));
temporal_pscs_subjavg = squeeze(mean(temporal_pscs, 1, 'omitnan'));

visual_pscs_subjhemiavg = squeeze(mean(visual_pscs_subjavg, 1, 'omitnan'));
auditory_pscs_subjhemiavg = squeeze(mean(auditory_pscs_subjavg, 1, 'omitnan'));
spatial_pscs_subjhemiavg = squeeze(mean(spatial_pscs_subjavg, 1, 'omitnan'));
temporal_pscs_subjhemiavg = squeeze(mean(temporal_pscs_subjavg, 1, 'omitnan'));

%% Plot PSC per ROI/hemi
set(groot, 'DefaultAxesTickLabelInterpreter', 'none');

% subj, hemisphere, ROI PSCs
lh_ROIs = cellfun(@(x) [x '_lh'], ROIs, 'UniformOutput', false);
rh_ROIs = cellfun(@(x) [x '_rh'], ROIs, 'UniformOutput', false);
hemi_ROIs = [lh_ROIs; rh_ROIs];
hemi_ROIs = hemi_ROIs(:)'; % interleave lh and rh
max_psc = max(psc_results,[],'all');
min_psc = min(psc_results,[],'all');

for ii = 1:N_contrasts

    psc_onecontrast = squeeze(psc_results(:,:,ii,:));
    psc_flat = psc_onecontrast(:);
    xs = repelem(1:(N_ROIs*2), N);
    xlabs = hemi_ROIs;

    figure;
    swarmchart(xs, psc_flat, 25, 'blue', 'filled');
    xticks(1:(N_ROIs*2));
    xticklabels(xlabs);
    yline(0);
    xlabel('ROIs');
    ylabel('% Signal Change');
    title(contrasts{ii});

    % Plot means too
    hold on;
    psc_subjmean_onecontrast = squeeze(psc_subjavg(:,ii,:));
    pscs_subjmean_flat = psc_subjmean_onecontrast(:);
    xs_avg = 1:(N_ROIs*2);
    swarmchart(xs_avg, pscs_subjmean_flat, 50, 'red', 'filled');

    ylim([min_psc, max_psc]);

end



%% Plot ROIs on spatial-temporal and visual-auditory axes
spatial_temporal = spatial_pscs - temporal_pscs;
visual_auditory = visual_pscs - auditory_pscs;
ROIs_per_modality = [4,4,6]; % MD, visual, auditory
modality_colors = {'g', 'b', 'r'};
modality_shapes = {'^', 'o', 'd'};

figure;
yline(0);
xline(0);
hold on;
handles = {};
ROI_count = 1;
for mm = 1:length(ROIs_per_modality)
    ROI_inds = ROI_count:ROI_count+ROIs_per_modality(mm)-1;
    modality_spatial_temporal = squeeze(spatial_temporal(:,:,ROI_inds));
    mean_modality_st = mean(modality_spatial_temporal, 1, 'omitnan');
    modality_visual_auditory = squeeze(visual_auditory(:,:,ROI_inds));
    mean_modality_va = mean(modality_visual_auditory, 1, 'omitnan');
    coords = [mean_modality_va(:), mean_modality_st(:)];
    handles{mm} = scatter(coords(:,1), coords(:,2), 100, modality_colors{mm}, modality_shapes{mm});
    hold on;
    ROI_count = ROI_count + ROIs_per_modality(mm);
    disptable = table(repelem(ROIs(ROI_inds),2)', coords(:,1), coords(:,2), 'VariableNames', {'ROI', 'v-a', 's-t'});
    disp(disptable);
end
xlabel('Visual-Auditory (PSC)');
ylabel('Spatial-Temporal (PSC)');
legend([handles{1}, handles{2}, handles{3}], 'Multiple Demand', 'Visual', 'Auditory')


%% Plot ROIs on passive visual-auditory axes
visual_auditory_passive = squeeze(psc_results(:,:,5,:)) - squeeze(psc_results(:,:,6,:));
ROIs_per_modality = [4,4,6]; % MD, visual, auditory
modality_colors = {'g', 'b', 'r'};
modality_shapes = {'^', 'o', 'd'};

figure;
yline(0);
xline(0);
hold on;
handles = {};
ROI_count = 1;
for mm = 1:length(ROIs_per_modality)
    ROI_inds = ROI_count:ROI_count+ROIs_per_modality(mm)-1;
    modality_visual_auditory_passive = visual_auditory_passive(:,:,ROI_inds);
    mean_modality_va_passive = squeeze(mean(modality_visual_auditory_passive, 1, 'omitnan'));
    handles{mm} = scatter(mean_modality_va_passive(:), zeros(length(mean_modality_va_passive(:)),1), 100, modality_colors{mm}, modality_shapes{mm});
    hold on;
    ROI_count = ROI_count + ROIs_per_modality(mm);
    disptable = table(repelem(ROIs(ROI_inds),2)', mean_modality_va_passive(:), 'VariableNames', {'ROI', 'pV-pA'});
    disp(disptable);
end
xlabel('Visual-Auditory (PSC)');
legend([handles{1}, handles{2}, handles{3}], 'Multiple Demand', 'Visual', 'Auditory')




