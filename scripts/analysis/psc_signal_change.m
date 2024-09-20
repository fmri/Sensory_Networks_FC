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
annot_file_suffix = '_ROIs_pVisMod.annot';

contrasts = {'sV-pV', 'tV-pV', 'sA-pA', 'tA-pA'};
reject_contrast = {{}, {'RT', 'LA'}, {'LN', 'GG', 'TP'}, {}}; % These subjs had below 55% accuracy on these tasks and should be rejected from analysis in those conditions
N_contrasts = length(contrasts);

ROIs = {'aINS', 'preSMA', 'ppreCun', 'dACC', ... % multisensory
    'sPCS', 'iPCS', 'midIFS', 'pVis_mod' ... % visual
    'tgPCS', 'cIFSG', 'pAud', 'CO', 'FO', 'cmSFG'}; % auditory
N_ROIs = length(ROIs);
load('/projectnb/somerslab/tom/projects/spacetime_network/data/missing_ROIs.mat', 'missing_ROIs');

hemis = {'lh', 'rh'};

psc_results = nan(N, 2, N_contrasts, N_ROIs); % subjs x hemis x contrasts x ROIs = 20 x 2 x 4 x 14

%% Loop through subj/hemisphere/contrast/ROI and calculate % signal change in ROI

for ss = 1:N % for subj

    subjCode = subjCodes{ss};

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
                psc_results(ss, hh, cc, :) = nan;
                continue;
            end

            % load subj contrast percent signal change (psc) data
            contrast_fpath = [subjdata_dir subjCode '/bold/spacetime_contrasts_' hemi '_newcondfiles/' contrast '/cespct.nii.gz'];
            psc_data = MRIread(contrast_fpath);

            for rr = 1:N_ROIs % for ROI

                ROI = ROIs{rr};

                if ismember([subjCode '_' ROI '_' hemi], missing_ROIs) 
                    psc_results(ss, hh, cc, rr) = nan;
                    continue;
                end

                ROI_index = find(ismember(annot_ctable.struct_names, ROI));
                ROI_label = annot_ctable.table(ROI_index,5);
                ROI_mask = annot_labels == ROI_label;
                assert(sum(ROI_mask)>=100, ['Subj ' subjCode ' ' hemi ' ' ROI ' is below 100 vertices']);
                
                % Get mean psc data from ROI mask
                psc_results(ss,hh,cc,rr) = mean(psc_data.vol(ROI_mask)); %%%% MAKE SURE THESE INDS LINE UP CORRECTLY
                assert(~isnan(psc_results(ss,hh,cc,rr)), ['PSC for subj ' subjCode ' ' hemi ' ' ROI ' computed as NaN'])

            end % end ROI

        end % end contrast
    end % end hemisphere

end % end subj


%% Average results over subjs
psc_subjavg = squeeze(mean(psc_results, 1, 'omitnan'));
psc_subjhemiavg = squeeze(mean(psc_results, [1,2], 'omitnan'));

%% Plot
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










