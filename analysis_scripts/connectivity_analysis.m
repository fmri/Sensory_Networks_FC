%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract the ROI timeseries data for each
% subject and use correlation analysis to calculate functional connectivity
% as well as hierarchical clustering to do network analysis
% Tom Possidente - June 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccc;

%% Setup analysis parameters
resting_state = false;
hierarchical_clustering = true;
plot_individual_connmats = false;
save_out = true;

if resting_state
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_resting_state/results/preprocessing/';
    N = 13;
    subjIDs = {'RR', 'MK', 'LA', 'TP', 'NM', 'SL', 'UV', 'PQ', 'KQ', 'LN', 'PT', 'PL', 'NS'};
    Ncond = 1; % one condition = resting state
    conditions = {'rest'};
    ROI_data = cell(N,1);
    %ROI_trial_inds = cell(N,1);
    fourth_ROI_check = 'all_ROIs.CO (L)';
    extra_ROIs_end = 4;
else
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_task_nofmap/results/preprocessing/';
    %ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_all_preproc_alt_fmap_order/results/preprocessing/';
    N = 16;
    subjIDs = {'RR', 'MM', 'PP', 'MK', 'LA', 'TP', 'NM', 'SL', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS'};
    Ncond = 10;
    conditions = {'Fixation' , 'Passive_Auditory', 'Passive_Tactile', 'Passive_Visual', 'Spatial_Auditory',...
        'Spatial_Tactile', 'Spatial_Visual', 'Temporal_Auditory', 'Temporal_Tactile', 'Temporal_Visual'};
    ROI_data = cell(N,Ncond);
    ROI_trial_inds = cell(N,1);
    fourth_ROI_check = 'all_ROIs.CO (L)'; %%%%
    extra_ROIs_end = 13; %%%%
end

%% Setup ROIs

aud_ROIs_use = {'tgPCS', 'cIFS/G', 'FO', 'CO', 'pAud'};
vis_ROIs_use = {'pVis', 'preSMA-V', 'SPCS', 'IPCS', 'midIFS'};
mult_ROIs_use = {'cmSFG_mult', 'Ins_mult'};

desired_order = {'tgPCS (L)', 'FO (L)', 'CO (L)', 'cIFS_G (L)', 'pAud (L)', 'tgPCS (R)', 'FO (R)', 'CO (R)', ...
    'cIFS_G (R)', 'pAud (R)',...
    'SPCS (L)', 'IPCS (L)', 'midIFS (L)', 'preSMA-V (L)', 'pVis (L)', 'SPCS (R)', 'IPCS (R)', 'midIFS (R)', ...
    'preSMA-V (R)', 'pVis (R)',...
    'cmSFG_mult (L)', 'Ins_mult (L)', 'cmSFG_mult (R)', 'Ins_mult (R)'};

% desired_order = {'tgPCS (L)', 'cIFS_G (L)', 'pAud (L)', 'tgPCS (R)', 'cIFS_G (R)', 'pAud (R)',...
%     'SPCS (L)', 'IPCS (L)', 'pVis (L)', 'SPCS (R)', 'IPCS (R)', 'pVis (R)',...
%     'cmSFG_mult (L)', 'Ins_mult (L)', 'cmSFG_mult (R)', 'Ins_mult (R)'};

%% Get missing ROI data
load('/projectnb/somerslab/tom/projects/spacetime_network/data/missing_ROIs.mat', 'missing_ROIs_allsubj', 'subjCodes');

%% Get ROI data

files = {dir(ROI_dataDir).name};

for ss = 1:N
    subjnum = sprintf( '%03d', ss ) ;
    subjfiles = files(contains(files, ['ROI_Subject' subjnum '_Condition']) & ~contains(files, '000'));

    condcount = 0;
    for ff = 1:length(subjfiles)

        load([ROI_dataDir subjfiles{ff}], 'names', 'data', 'conditionname', 'conditionweights');
        %data2 = data;
        if resting_state
            assert(strcmp(conditionname,'rest'), ['Condition name is not "rest" for subj ' num2str(ss)])
            condcount = condcount + 1;
            cond_ind = 1;
        else
            if ~ismember(conditionname, conditions)
                disp(['Subj ' num2str(ss) ' has unrecognized condition name ' conditionname ' ...skipping...'])
                continue
            else
                condcount = condcount + 1;
                [~,cond_ind] = ismember(conditionname, conditions);
            end

            % Index for condition data only (boxcar filter)
            conditionweight = conditionweights{2};
            conditionweight(conditionweight>0) = 1;
            data = cellfun(@(x) x.*conditionweight, data, 'UniformOutput',false); % multiply condition on/off by ROI data to get condition data only
            %trial_inds_pre = cellfun(@(x) find(abs(x)>0), data, 'UniformOutput',false);
            %ROI_trial_inds{ss,cond_ind} = cellfun(@(x) diff(x)>1, trial_inds_pre, 'UniformOutput',false);
            data = cellfun(@(x) x(abs(x)>0), data, 'UniformOutput',false);

            % Alternative indexing for condition data only (gaussian instead of boxcar)
            %data2 = cellfun(@(x) x.*conditionweights{2}, data2, 'UniformOutput',false); % multiply condition on/off by ROI data to get condition data only
            %data2 = cellfun(@(x) x(abs(x)>0), data2, 'UniformOutput',false);

        end

        assert( strcmp(names{3},'CSF') & strcmp(names{4}, fourth_ROI_check) , ['Unexpected ROI order for subj: ' num2str(ff)]);
        data = data(4:end-extra_ROIs_end);
        %ROI_trial_inds{ss,cond_ind} = ROI_trial_inds{ss,cond_ind}(4:end-extra_ROIs_end);
        %data2 = data2(4:end-extra_ROIs_end);
        names = names(4:end-extra_ROIs_end);
        names_clean = cellfun(@(f) f(10:end), names, 'UniformOutput', false);

        % Reorder ROIs
        [ROIs_match, reorder_inds] = ismember(desired_order, names_clean);
        %assert(sum(ROIs_match)==24, ['incorrect ROIs for subj ' num2str(ss) ' condition ' conditions{cond_ind}]);
        names = names_clean(reorder_inds);
        data = data(reorder_inds);
        ROI_data{ss,cond_ind} = data;
        %ROI_data2{ss,cond_ind} = data2;

    end

    if resting_state
        assert(condcount==1, ['Unexpected number of conditions for subj ' num2str(ss)]);
    else
        assert(condcount==10, ['Unexpected number of conditions for subj ' num2str(ss)]);
    end

end


%% Connectivity correlations per subj

N_ROIs = length(names);
connmats = NaN(N_ROIs, N_ROIs, N, Ncond);
%connmats2 = NaN(N_ROIs, N_ROIs, N, Ncond);
pvals = NaN(N_ROIs, N_ROIs, N, Ncond);

for ss = 1:N
    subjID = subjIDs{ss};
    for cc = 1:Ncond
        data = ROI_data{ss,cc};
        %data2 = ROI_data2{ss,cc};
        for rr1 = 1:N_ROIs
            ROI1 = data{rr1};
            %ROI1_2 = data2{rr1};
            for rr2 = 1:N_ROIs
                ROI2 = data{rr2};
                %ROI2_2 = data2{rr2};
                %if resting_state
                missing_ROIs = missing_ROIs_allsubj{ismember(subjIDs,subjID)};
                if isempty(missing_ROIs)
                    [connmats(rr1,rr2,ss,cc), pvals(rr1,rr2,ss,cc)] = corr(ROI1, ROI2);
                elseif any(ismember({desired_order{rr1}, desired_order{rr2}}, missing_ROIs)) % if subj is missing ROI, connmat is nan
                    connmats(rr1,rr2,ss,cc) = NaN;
                    pvals(rr1,rr2,ss,cc) = NaN;
                    disp(['Subj ' num2str(ss) ' (' subjID ') missing ROIs: ' missing_ROIs]);
                else
                    [connmats(rr1,rr2,ss,cc), pvals(rr1,rr2,ss,cc)] = corr(ROI1, ROI2);
                end
                 
                    %[connmats2(rr1,rr2,ss,cc), ~] = corr(ROI1_2, ROI2_2);
                %else
                %    ROI1_inds = [find(ROI_trial_inds{ss,cc}{rr1}); length(ROI1)]; %% This code separates trials, correlates trial signals then averages the correlation coeffs
                %    ROI2_inds = [find(ROI_trial_inds{ss,cc}{rr2}); length(ROI1)];
                %    N_trials = length(find(ROI1_inds));
                %    trial_corrs = NaN(N_trials,1);
                %    trial_pvals = NaN(N_trials,1);
                %    for tt = 1:N_trials
                %        trial_data1 = ROI1(ROI1_inds(tt)-ROI1_inds(1)+1:ROI1_inds(tt));
                %        trial_data2 = ROI2(ROI2_inds(tt)-ROI2_inds(1)+1:ROI2_inds(tt));
                %        [trial_corrs(tt), trial_pvals(tt)] = corr(trial_data1, trial_data2);
                %    end
                %    connmats(rr1,rr2,ss,cc) = mean(trial_corrs);
                %    pvals(rr1,rr2,ss) = max(trial_pvals);
                %end
            end
        end
    end
end

%% Connectivity group
connmat_group = mean(connmats, 3, 'omitnan');
%connmat_group2 = mean(connmats2,3);

%% Calculate hierarchical clustering
if hierarchical_clustering
    for cc = 1:Ncond

        distance_measure = 1-abs(connmat_group(:,:,1,cc));
        linkage_cluster = linkage(distance_measure);

        figure;
        dendrogram(linkage_cluster,'Labels',names);
        title(['Hierarchical Clustering Condition ' conditions{cc}]);

    end
end

%% Plot group level connectivity matrices

connmat_group_nandiag = NaN(size(connmat_group));
for mm = 1:size(connmat_group,4)
    connmat_group_nandiag(:,:,1,mm) = connmat_group(:,:,1,mm) - diag(diag(connmat_group(:,:,1,mm))) + diag(NaN(size(connmat_group,1),1));
end

for cc = 1:Ncond
    figure;
    heatmap(names, names, connmat_group_nandiag(:,:,1,cc), 'Colormap', turbo, 'ColorLimits', [min(connmat_group_nandiag(),[], 'all'),max(connmat_group_nandiag,[],'all')], 'FontSize',16)
    title(['Group Average Conn Mat | Condition ' conditions{cc}]);
end

% for cc = 1:Ncond
%     figure;
%
%     %[min(connmat_group,[], 'all'),max(connmat_group,[],'all')]
%     heatmap(names, names, connmat_group2(:,:,1,cc), 'Colormap', turbo, 'ColorLimits', [-1,1], 'FontSize',16)
%     title(['Group Average Conn Mat 2 | Condition ' conditions{cc}])
% end

%% Plot group level connectivity difference matrices 
if ~resting_state
    ztrans_connmat_group = atanh(connmat_group_nandiag);% fisher z-transoform the rs
    diffs1 = [2,4, 5,8,  7,10, 5,7];
    diffs2 = [1,1, 2,2 , 4,4,  8,10];
        conditions = {'Fixation' , 'Passive_Auditory', 'Passive_Tactile', 'Passive_Visual', 'Spatial_Auditory',...
            'Spatial_Tactile', 'Spatial_Visual', 'Temporal_Auditory', 'Temporal_Tactile', 'Temporal_Visual'};
    for cc = 1:length(diffs1)
        figure;
        diff = ztrans_connmat_group(:,:,1,diffs1(cc))-ztrans_connmat_group(:,:,1,diffs2(cc));
        heatmap(names, names,  diff, 'Colormap', turbo, 'ColorLimits', [min(diff(),[], 'all'),max(diff,[],'all')], 'FontSize',16)
        title([conditions{diffs1(cc)} '-' conditions{diffs2(cc)}])
    end
end

%% Plot individual connectivity matrices

%sig_inds = pvals >= 0.01; % blank out correlations that are not significant
%connmats(sig_inds) = nan;
connmats_nandiag = NaN(size(connmats));
for cc = 1:Ncond
    for mm = 1:size(connmats,3)
        connmats_nandiag(:,:,mm,cc) = connmats(:,:,mm,cc) - diag(diag(connmats(:,:,mm,cc))) + diag(NaN(size(connmats,1),1));
    end
end

%[min(connmats(:,:,pp),[], 'all'),max(connmats(:,:,pp),[],'all')]
if plot_individual_connmats
    for cc = 1:Ncond
        figure;
        count = 0;
        for pp1 = 1:4
            for pp2 = 1:4
                count = count + 1;
                subplot(4,4,count);
                if count <= N
                    heatmap(names, names, connmats_nandiag(:,:,count,cc), 'Colormap', turbo, 'ColorLimits', [min(connmats_nandiag(:,:,count,cc),[], 'all', 'omitnan'),max(connmats_nandiag(:,:,count,cc),[],'all','omitnan')])
                    title(['Conn Mat Subj ' num2str(count) ' condition ' conditions{cc}])
                end
            end
        end
    end
end


%% Save out connectivity results
if save_out

    if resting_state
        save('rs_group_connmat.mat', 'names', 'connmat_group_nandiag');
        save('rs_indiv_connmats.mat', 'names', 'connmats_nandiag');
    else
        save('task_group_connmat.mat', 'names', 'connmat_group_nandiag');
        save('task_indiv_connmats.mat', 'names', 'connmats_nandiag');
    end

end
