%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to quantify the overlap of the
%%% significant sensory drive and working memory areas for all ROIs using
%%% the localizer data from the spacetime experiment
%%%
%%% Tom Possidente - December 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load in tstat data
load('/projectnb/somerslab/tom/projects/spacetime_network/data/probROI_sensory_WM_sigstats_localizer.mat', ...
    'ROI_groupavg_tstats_MD', 'ROI_groupavg_tstats', 'RAS_coords_MD', 'RAS_coords', 'tstats_act', ...
    'tstats_pass', 'tstats_pass_MD', 'tstats_act_MD', 'ROI_list', 'subjCodes', 'hemis', ...
    'active_contrast_list', 'passive_contrast_list');
% I saved these with the wrong variable names, so switch to correct names
sigstats_WM = tstats_act;
sigstats_sensory = tstats_pass;
sigstats_WM_MD = tstats_act_MD;
sigstats_sensory_MD = tstats_pass_MD;

N_subjs = length(subjCodes);
N_ROIs = length(ROI_list);

%% Initialize Variables
load('/projectnb/somerslab/tom/projects/spacetime_network/data/missing_ROIs.mat', 'missing_ROIs');
load('/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/replacement_ROI_list.mat', 'replacement_ROIs');
overlap = nan(N_subjs, N_ROIs-3, 2, 3); % subjs x ROIs x hemis x [percent active, percent passive, percent overlap]
overlap_MD = nan(N_subjs, N_ROIs-10, 2, 3, 2); % multiple demand ROIs will have 2 overlaps, one for auditory, one for visual
sigthresh = 0.05;

%% Loop through subjs and find overlap % of significant WM and sensory regions within ROI
% %-log10(p)
count = 0;
bad_count = 0;

for rr = 1:N_ROIs
    for hh = 1:2

        if rr<=10
            ROI_RAS_coords = RAS_coords{rr,hh};
            N_subjs_ROI = size(sigstats_sensory{rr,hh},2);
        else
            ROI_RAS_coords = RAS_coords_MD{rr-10,hh};
            N_subjs_ROI = size(sigstats_sensory_MD{rr-10,hh},2);
        end

        for ss = 1:N_subjs_ROI
            ROI_code = [subjCodes{ss} '_' ROI_list{rr} '_' hemis{hh}];
            if ismember(ROI_code, replacement_ROIs) || ismember(ROI_code, missing_ROIs)
                continue;
            end
            if rr<=10
                count = count + 1;
                sensory_sig = sigstats_sensory{rr,hh}(:,ss) > -log10(sigthresh);
                WM_sig = sigstats_WM{rr,hh}(:,ss) > -log10(sigthresh);
                total_sig = sum(sensory_sig) + sum(WM_sig);

                if sum(sensory_sig)==0
                    bad_count = bad_count + 1;
                    disp(['subj ' num2str(ss) ' ' ROI_list{rr} ' ' hemis{hh} ' has no significant sensory vertices' ])
                    continue
                end

                if sum(WM_sig)==0
                    bad_count = bad_count + 1;
                    disp(['subj ' num2str(ss) ' ' ROI_list{rr} ' ' hemis{hh} ' has no significant WM vertices' ])
                    continue
                end

                overlap(ss,rr,hh,1) = sum(sensory_sig)*2 / total_sig;
                overlap(ss,rr,hh,2) = sum(sensory_sig & WM_sig)*2 / total_sig;
                overlap(ss,rr,hh,3) = sum(WM_sig)*2 / total_sig;
            else
                % Visual
                count = count + 1;
                sensory_sig = sigstats_sensory_MD{rr-10 ,hh, 1}(:,ss) > -log10(sigthresh);
                WM_sig = sigstats_WM_MD{rr-10,hh, 1}(:,ss) > -log10(sigthresh);
                total_sig = sum(sensory_sig) + sum(WM_sig);

                if sum(sensory_sig)==0
                    bad_count = bad_count + 1;
                    disp(['subj ' num2str(ss) ' ' ROI_list{rr} ' ' hemis{hh} ' has no significant sensory vertices' ])
                    continue
                end

                if sum(WM_sig)==0
                    bad_count = bad_count + 1;
                    disp(['subj ' num2str(ss) ' ' ROI_list{rr} ' ' hemis{hh} ' has no significant WM vertices' ])
                    continue
                end

                overlap_MD(ss,rr-10,hh,1,1) = sum(sensory_sig)*2 / total_sig;
                overlap_MD(ss,rr-10,hh,2,1) = sum(sensory_sig & WM_sig)*2 / total_sig;
                overlap_MD(ss,rr-10,hh,3,1) = sum(WM_sig)*2 / total_sig;

                % Auditory
                count = count + 1;
                sensory_sig = sigstats_sensory_MD{rr-10 ,hh, 2}(:,ss) > -log10(sigthresh);
                WM_sig = sigstats_WM_MD{rr-10,hh, 2}(:,ss) > -log10(sigthresh);
                total_sig = sum(sensory_sig) + sum(WM_sig);

                if sum(sensory_sig)==0
                    bad_count = bad_count + 1;
                    disp(['subj ' num2str(ss) ' ' ROI_list{rr} ' ' hemis{hh} ' has no significant sensory vertices' ])
                    continue
                end

                if sum(WM_sig)==0
                    bad_count = bad_count + 1;
                    disp(['subj ' num2str(ss) ' ' ROI_list{rr} ' ' hemis{hh} ' has no significant WM vertices' ])
                    continue
                end

                overlap_MD(ss,rr-10,hh,1,2) = sum(sensory_sig)*2 / total_sig;
                overlap_MD(ss,rr-10,hh,2,2) = sum(sensory_sig & WM_sig)*2 / total_sig;
                overlap_MD(ss,rr-10,hh,3,2) = sum(WM_sig)*2 / total_sig;

            end

        end
    end
end


%% Plot

for rr = 1:N_ROIs

    if rr<=10
        sensory_perc = overlap(:,rr,:,1);
        overlap_perc = overlap(:,rr,:,2);
        WM_perc = overlap(:,rr,:,3);
        
        disp([ROI_list{rr} ' DICE: ' num2str(mean(overlap_perc(:),1,'omitnan'))] );

        figure;
        histogram(overlap_perc(:),5);
        title([ROI_list{rr} ' overlap']);
        xlim([0,1]);

    else
        % Visual
        sensory_perc = overlap_MD(:,rr-10,:,1,1);
        overlap_perc = overlap_MD(:,rr-10,:,2,1);
        WM_perc = overlap_MD(:,rr-10,:,3,1);

        figure;
        histogram(overlap_perc(:),5);
        title([ROI_list{rr} ' visual overlap']);
        xlim([0,1]);

        disp([ROI_list{rr} ' visual DICE: ' num2str(mean(overlap_perc(:),1,'omitnan'))] );

        % Auditory
        sensory_perc = overlap_MD(:,rr-10,:,1,2);
        overlap_perc = overlap_MD(:,rr-10,:,2,2);
        WM_perc = overlap_MD(:,rr-10,:,3,2);

        figure;
        histogram(overlap_perc(:),5);
        title([ROI_list{rr} ' auditory overlap']);
        xlim([0,1]);

        disp([ROI_list{rr} ' auditory DICE: ' num2str(mean(overlap_perc(:),1,'omitnan'))] );

    end

end






