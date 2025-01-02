%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to find the weighted centroids of the
%%% sensory drive tstat maps and working memory tstat maps for all ROIs using
%%% the localizer data from the spacetime experiment
%%%
%%% Tom Possidente - December 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Load in tstat data
load('/projectnb/somerslab/tom/projects/spacetime_network/data/probROI_sensory_WM_tstats_localizer.mat', ...
    'ROI_groupavg_tstats_MD', 'ROI_groupavg_tstats', 'RAS_coords_MD', 'RAS_coords', 'tstats_act', ...
    'tstats_pass', 'tstats_pass_MD', 'tstats_act_MD', 'ROI_list', 'subjCodes', 'hemis', ...
    'active_contrast_list', 'passive_contrast_list')
N_subjs = length(subjCodes);
N_ROIs = length(ROI_list);

%% Initialize Variables
load('/projectnb/somerslab/tom/projects/spacetime_network/data/missing_ROIs.mat', 'missing_ROIs');
load('/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/replacement_ROI_list.mat', 'replacement_ROIs');
wc_RAS_sensory = nan(N_subjs, N_ROIs-3, 2, 3); % subjs x ROIs x hemis x RAS coordinates
wc_RAS_WM = nan(N_subjs, N_ROIs-3, 2, 3); % wc = weighted centroid, RAS = right, anterior, superior coordinate system
wc_RAS_sensory_MD = nan(N_subjs, N_ROIs-10, 2, 3, 2); % multiple demand ROIs will have 2 sets of centroids, one for auditory, one for visual
wc_RAS_WM_MD = nan(N_subjs, N_ROIs-10, 2, 3, 2);
tthresh = 0;
mm_thresh = 0;

%% Loop through ROIs/subjs/hemispheres and calculate weighted centroid for both tstat maps (sensory & WM)

for rr = 1:N_ROIs
    for hh = 1:2

        if rr<=10
            ROI_RAS_coords = RAS_coords{rr,hh};
            N_subjs_ROI = size(tstats_pass{rr,hh},2);
        else
            ROI_RAS_coords = RAS_coords_MD{rr-10,hh};
            N_subjs_ROI = size(tstats_pass_MD{rr-10,hh},2);
        end

        for ss = 1:N_subjs_ROI
            ROI_code = [subjCodes{ss} '_' ROI_list{rr} '_' hemis{hh}];
            if ismember(ROI_code, replacement_ROIs) || ismember(ROI_code, missing_ROIs)
                continue;
            end
            if rr<=10 % sensory-biased ROIs
                sensory_tstats = tstats_pass{rr,hh}(:,ss);
                sensory_tstats(sensory_tstats<tthresh) = 0;
                weight_sum = sum(sensory_tstats);
                wc_RAS_sensory(ss,rr,hh,1) = sum(sensory_tstats.*ROI_RAS_coords{:,2}) / weight_sum;
                wc_RAS_sensory(ss,rr,hh,2) = sum(sensory_tstats.*ROI_RAS_coords{:,3}) / weight_sum;
                wc_RAS_sensory(ss,rr,hh,3) = sum(sensory_tstats.*ROI_RAS_coords{:,4}) / weight_sum;

                WM_tstats = tstats_act{rr,hh}(:,ss);
                WM_tstats(WM_tstats<tthresh) = 0;
                weight_sum = sum(WM_tstats);
                wc_RAS_WM(ss,rr,hh,1) = sum(WM_tstats.*ROI_RAS_coords{:,2}) / weight_sum;
                wc_RAS_WM(ss,rr,hh,2) = sum(WM_tstats.*ROI_RAS_coords{:,3}) / weight_sum;
                wc_RAS_WM(ss,rr,hh,3) = sum(WM_tstats.*ROI_RAS_coords{:,4}) / weight_sum;
            else % multidemand ROIs
                sensory_tstats_vis = tstats_pass_MD{rr-10,hh,1}(:,ss);
                sensory_tstats_vis(sensory_tstats_vis<tthresh) = 0;
                weight_sum = sum(sensory_tstats_vis);
                wc_RAS_sensory_MD(ss,rr-10,hh,1,1) = sum(sensory_tstats_vis.*ROI_RAS_coords{:,2}) / weight_sum;
                wc_RAS_sensory_MD(ss,rr-10,hh,2,1) = sum(sensory_tstats_vis.*ROI_RAS_coords{:,3}) / weight_sum;
                wc_RAS_sensory_MD(ss,rr-10,hh,3,1) = sum(sensory_tstats_vis.*ROI_RAS_coords{:,4}) / weight_sum;

                sensory_tstats_aud = tstats_pass_MD{rr-10,hh,2}(:,ss);
                sensory_tstats_aud(sensory_tstats_aud<tthresh) = 0;
                weight_sum = sum(sensory_tstats_aud);
                wc_RAS_sensory_MD(ss,rr-10,hh,1,2) = sum(sensory_tstats_aud.*ROI_RAS_coords{:,2}) / weight_sum;
                wc_RAS_sensory_MD(ss,rr-10,hh,2,2) = sum(sensory_tstats_aud.*ROI_RAS_coords{:,3}) / weight_sum;
                wc_RAS_sensory_MD(ss,rr-10,hh,3,2) = sum(sensory_tstats_aud.*ROI_RAS_coords{:,4}) / weight_sum;

                WM_tstats_vis = tstats_act_MD{rr-10,hh,1}(:,ss);
                WM_tstats_vis(WM_tstats_vis<tthresh) = 0;
                weight_sum = sum(WM_tstats_vis);
                wc_RAS_WM_MD(ss,rr-10,hh,1,1) = sum(WM_tstats_vis.*ROI_RAS_coords{:,2}) / weight_sum;
                wc_RAS_WM_MD(ss,rr-10,hh,2,1) = sum(WM_tstats_vis.*ROI_RAS_coords{:,3}) / weight_sum;
                wc_RAS_WM_MD(ss,rr-10,hh,3,1) = sum(WM_tstats_vis.*ROI_RAS_coords{:,4}) / weight_sum;

                WM_tstats_aud = tstats_act_MD{rr-10,hh,2}(:,ss);
                WM_tstats_aud(WM_tstats_aud<tthresh) = 0;
                weight_sum = sum(WM_tstats_aud);
                wc_RAS_WM_MD(ss,rr-10,hh,1,2) = sum(WM_tstats_aud.*ROI_RAS_coords{:,2}) / weight_sum;
                wc_RAS_WM_MD(ss,rr-10,hh,2,2) = sum(WM_tstats_aud.*ROI_RAS_coords{:,3}) / weight_sum;
                wc_RAS_WM_MD(ss,rr-10,hh,3,2) = sum(WM_tstats_aud.*ROI_RAS_coords{:,4}) / weight_sum;
            end
        end
    end
end


%% Plot results

for rr = 1:N_ROIs

    % Anterior-Posterior, Inferior-Superior coordinate analysis
    if rr <= 10
        wc_WM_ant = squeeze(wc_RAS_WM(:,rr,:,2));
        wc_sensory_ant = squeeze(wc_RAS_sensory(:,rr,:,2));
        wc_diff_ant = wc_WM_ant(:) - wc_sensory_ant(:);
        wc_WM_sup = squeeze(wc_RAS_WM(:,rr,:,3));
        wc_sensory_sup = squeeze(wc_RAS_sensory(:,rr,:,3));
        wc_diff_sup = wc_WM_sup(:) - wc_sensory_sup(:);
        wc_WM_right = squeeze(wc_RAS_WM(:,rr,:,1));
        wc_sensory_right = squeeze(wc_RAS_sensory(:,rr,:,1));
        wc_diff_right = wc_WM_right(:) - wc_sensory_right(:);
        wc_dist = sqrt(  (wc_diff_ant.^2 + wc_diff_sup.^2 + wc_diff_right.^2) );
        ang = atan(wc_diff_sup./wc_diff_ant);
        ang(wc_diff_ant<0) = ang(wc_diff_ant<0) - pi; % if anterior coordinate is negative, add 180 to angle to put it in the correct quadrant
        [h,p,ci,stats] = ttest(wc_diff_ant(:), mm_thresh, 'tail', 'right');
        disp([ROI_list{rr} ' p=' num2str(round(p,5))]);
        disp(['Mean RAS diffs: ' num2str(nanmean(wc_diff_right)) ', ' num2str(nanmean(wc_diff_ant)) ', ' num2str(nanmean(wc_diff_sup)) ] );
        
        wc_diff_sup_mean = mean(wc_diff_sup,1,'omitnan');
        wc_diff_ant_mean = mean(wc_diff_ant, 1, 'omitnan');
        wc_diff_right_mean = mean(wc_diff_right, 1, 'omitnan');
        wc_dist_mean = sqrt(  (wc_diff_ant_mean.^2 + wc_diff_sup_mean.^2 + wc_diff_right_mean.^2) );
        ang_mean = atan(wc_diff_sup_mean./wc_diff_ant_mean);
        if wc_diff_ant_mean<0; ang_mean = ang_mean - pi; end

        figure;
        for ss = 1:N_subjs*2
            if ss<=10
                polarplot([0 ang(ss)], [0 wc_dist(ss)], 'Marker', '.', 'Color', 'b');
            else
                polarplot([0 ang(ss)], [0 wc_dist(ss)], 'Marker', '.', 'Color', 'r');
            end
            hold on;
        end
        polarplot([0 ang_mean], [0 wc_dist_mean], 'Marker', '.', 'Color', 'g', 'LineWidth', 3);
        thetaticklabels({'Anterior', '', '', 'Superior', '', '', 'Posterior', '', '', 'Inferior', '', ''})
        title([ROI_list{rr} ' | p=' num2str(round(p,5))]);

        % figure;
        % scatter3()

    else
        wc_WM_ant = squeeze(wc_RAS_WM_MD(:,rr-10,:,2,1));
        wc_sensory_ant = squeeze(wc_RAS_sensory_MD(:,rr-10,:,2,1));
        wc_diff_ant = wc_WM_ant(:) - wc_sensory_ant(:);
        wc_WM_sup = squeeze(wc_RAS_WM_MD(:,rr-10,:,3,1));
        wc_sensory_sup = squeeze(wc_RAS_sensory_MD(:,rr-10,:,3,1));
        wc_diff_sup = wc_WM_sup(:) - wc_sensory_sup(:);
        wc_WM_right = squeeze(wc_RAS_WM_MD(:,rr-10,:,1,1));
        wc_sensory_right = squeeze(wc_RAS_sensory_MD(:,rr-10,:,1,1));
        wc_diff_right = wc_WM_right(:) - wc_sensory_right(:);
        wc_dist = sqrt(  (wc_diff_ant.^2 + wc_diff_sup.^2 + wc_diff_right.^2) );
        ang = atan(wc_diff_sup./wc_diff_ant);
        ang(wc_diff_ant<0) = ang(wc_diff_ant<0) - pi; % if anterior coordinate is negative, add 180 to angle to put it in the correct quadrant
        [h,p,ci,stats] = ttest(wc_diff_ant(:), mm_thresh, 'tail', 'right');
        disp([ROI_list{rr} ' visual p=' num2str(round(p,5))]);
        disp(['Mean RAS diffs: ' num2str(nanmean(wc_diff_right)) ', ' num2str(nanmean(wc_diff_ant)) ', ' num2str(nanmean(wc_diff_sup)) ] );

        wc_diff_sup_mean = mean(wc_diff_sup,1,'omitnan');
        wc_diff_ant_mean = mean(wc_diff_ant, 1, 'omitnan');
        wc_diff_right_mean = mean(wc_diff_right, 1, 'omitnan');
        wc_dist_mean = sqrt(  (wc_diff_ant_mean.^2 + wc_diff_sup_mean.^2 + wc_diff_right_mean.^2) );
        ang_mean = atan(wc_diff_sup_mean./wc_diff_ant_mean);
        if wc_diff_ant_mean<0; ang_mean = ang_mean - pi; end

        figure;
        for ss = 1:N_subjs*2
            if ss<=10
                polarplot([0 ang(ss)], [0 wc_dist(ss)], 'Marker', '.', 'Color', 'b');
            else
                polarplot([0 ang(ss)], [0 wc_dist(ss)], 'Marker', '.', 'Color', 'r');
            end            
            hold on;
        end
        polarplot([0 ang_mean], [0 wc_dist_mean], 'Marker', '.', 'Color', 'g', 'LineWidth', 3);
        thetaticklabels({'Anterior', '', '', 'Superior', '', '', 'Posterior', '', '', 'Inferior', '', ''})
        title([ROI_list{rr} ' visual | p=' num2str(round(p,5))]);

        wc_WM_ant = squeeze(wc_RAS_WM_MD(:,rr-10,:,2,2));
        wc_sensory_ant = squeeze(wc_RAS_sensory_MD(:,rr-10,:,2,2));
        wc_diff_ant = wc_WM_ant(:) - wc_sensory_ant(:);
        wc_WM_sup = squeeze(wc_RAS_WM_MD(:,rr-10,:,3,2));
        wc_sensory_sup = squeeze(wc_RAS_sensory_MD(:,rr-10,:,3,2));
        wc_diff_sup = wc_WM_sup(:) - wc_sensory_sup(:);
        wc_WM_right = squeeze(wc_RAS_WM_MD(:,rr-10,:,1,2));
        wc_sensory_right = squeeze(wc_RAS_sensory_MD(:,rr-10,:,1,2));
        wc_diff_right = wc_WM_right(:) - wc_sensory_right(:);
        wc_dist = sqrt(  (wc_diff_ant.^2 + wc_diff_sup.^2) );
        ang = atan(wc_diff_sup./wc_diff_ant);
        ang(wc_diff_ant<0) = ang(wc_diff_ant<0) - pi; % if anterior coordinate is negative, add 180 to angle to put it in the correct quadrant
        [h,p,ci,stats] = ttest(wc_diff_ant(:), mm_thresh, 'tail', 'right');
        disp([ROI_list{rr} ' visual p=' num2str(round(p,5))]);
        disp(['Mean RAS diffs: ' num2str(nanmean(wc_diff_right)) ', ' num2str(nanmean(wc_diff_ant)) ', ' num2str(nanmean(wc_diff_sup)) ] );

        wc_diff_sup_mean = mean(wc_diff_sup,1,'omitnan');
        wc_diff_ant_mean = mean(wc_diff_ant, 1, 'omitnan');
        wc_diff_right_mean = mean(wc_diff_right, 1, 'omitnan');
        wc_dist_mean = sqrt(  (wc_diff_ant_mean.^2 + wc_diff_sup_mean.^2 + wc_diff_right_mean.^2) );
        ang_mean = atan(wc_diff_sup_mean./wc_diff_ant_mean);
        if wc_diff_ant_mean<0; ang_mean = ang_mean - pi; end

        figure;
        for ss = 1:N_subjs*2
            if ss<=10
                polarplot([0 ang(ss)], [0 wc_dist(ss)], 'Marker', '.', 'Color', 'b');
            else
                polarplot([0 ang(ss)], [0 wc_dist(ss)], 'Marker', '.', 'Color', 'r');
            end
            hold on;
        end
        polarplot([0 ang_mean], [0 wc_dist_mean], 'Marker', '.', 'Color', 'g', 'LineWidth', 3);
        thetaticklabels({'Anterior', '', '', 'Superior', '', '', 'Posterior', '', '', 'Inferior', '', ''})
        title([ROI_list{rr} ' auditory | p=' num2str(round(p,5))]);

    end

    keyboard
end
