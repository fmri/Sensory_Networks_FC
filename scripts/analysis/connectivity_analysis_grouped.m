%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract the ROI timeseries data for each
% subject and use correlation analysis to calculate functional connectivity
% as well as hierarchical clustering to do network analysis
% Tom Possidente - June 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Setup analysis parameters
resting_state = false;
hierarchical_clustering = true;
plot_individual_connmats = false;
save_out = false;

if resting_state
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_resting_state/results/preprocessing/';
    subjCodes = {'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'AI', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'PT', 'PL', 'NS'};
    conditions = {'rest'};
    ROI_str = 'ROIs'; % "ROI_mod" for pVis without DO, "ROI" for pVis with DO
else % task data
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/preprocessing/';
    subjCodes = {'MM'	'PP' 'MK' 'AB' 'AD'	'LA' 'AE' 'TP' 'NM'	'AF' 'AG' 'GG' 'UV'	'PQ' 'KQ' 'LN' 'RT'	'PT' 'PL' 'NS'};
    conditions = {'f' , 'aP', 'tP', 'vP', 'aS', 'tS', 'vS', 'aT', 'tT', 'vT'};
    ROI_str = 'grouped_ROIs'; 
    reject_conditions = {{}, {}, {}, {}, {'LN', 'GG', 'TP'}, {}, {}, {'RT', 'LA'}, {}, {}}; % These subjs had below 55% accuracy on these tasks and should be rejected from analysis in those conditions
end

N = length(subjCodes);
Ncond = length(conditions);
ROI_data = cell(N,Ncond);
pVis_name = 'pVis'; % pVis_mod for pVis wihout DO, pVis for pVis with DO
first_ROI_check = [ROI_str '.MD_frontal (L)'];

%% Setup ROIs
aud_ROIs_use = {'abiased_frontal', 'pAud'};
vis_ROIs_use = {'vbiased_frontal', pVis_name};
mult_ROIs_use = {'MD_frontal'};

desired_order = {'abiased_frontal (L)' 'pAud (L)', 'abiased_frontal (R)', 'pAud (R)',...
    'vbiased_frontal (L)', [pVis_name ' (L)'], ...
    'vbiased_frontal (R)', [pVis_name ' (R)'], ...
    'MD_frontal (L)', 'MD_frontal (R)'};

%% Get ROI data

files = {dir(ROI_dataDir).name};

for ss = 1:N
    subjnum = sprintf( '%03d', ss ) ;
    subjfiles = files(contains(files, ['ROI_Subject' subjnum '_Condition']) & ~contains(files, '000'));
    
    condcount = 0;
    for ff = 1:length(subjfiles)

        load([ROI_dataDir subjfiles{ff}], 'names', 'data', 'conditionname', 'conditionweights');
        ROI_names_mask = cellfun(@(x) contains(x,ROI_str), names);
        names = names(ROI_names_mask);
        data = data(ROI_names_mask);

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
            conditionweight(conditionweight>0) = 1; % could get rid of this line to make it a tapered mask
            data = cellfun(@(x) x.*conditionweight, data, 'UniformOutput',false); % multiply condition on/off by ROI data to get condition data only
            data = cellfun(@(x) x(abs(x)>0), data, 'UniformOutput',false);

        end

        assert( strcmp(names{1}, first_ROI_check) , ['Unexpected ROI order for subj: ' num2str(ff)]);
        names_clean = cellfun(@(f) f(length(ROI_str)+2:end), names, 'UniformOutput', false);

        % Reorder ROIs
        [ROIs_match, reorder_inds] = ismember(desired_order, names_clean);
        names = names_clean(reorder_inds);
        data = data(reorder_inds);
        ROI_data{ss,cond_ind} = data;

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
pvals = NaN(N_ROIs, N_ROIs, N, Ncond);

for ss = 1:N
    subjID = subjCodes{ss};
    for cc = 1:Ncond
        if ismember(subjID, reject_conditions{cc})
            continue
        end
        data = ROI_data{ss,cc};
        for rr1 = 1:N_ROIs
            ROI1 = data{rr1};
            ROI1_name = replace([subjCodes{ss} '_' desired_order{rr1}], ' (L)', '_lh');
            ROI1_name = replace(ROI1_name, ' (R)', '_rh');
            for rr2 = 1:N_ROIs
                ROI2 = data{rr2};
                ROI2_name = replace([subjCodes{ss} '_' desired_order{rr2}], ' (L)', '_lh');
                ROI2_name = replace(ROI2_name, ' (R)', '_rh');
                [connmats(rr1,rr2,ss,cc), pvals(rr1,rr2,ss,cc)] = corr(ROI1, ROI2);
            end
        end
    end
end

%% Connectivity group
connmat_group = mean(connmats, 3, 'omitnan');

%% Calculate hierarchical clustering
if hierarchical_clustering
    for cc = 1:Ncond
        distance_measure = 1-abs(connmat_group(:,:,1,cc));
        linkage_cluster = linkage(distance_measure);

        figure;
        [h,t,outperm] = dendrogram(linkage_cluster,'Labels',names);
        title(['Hierarchical Clustering Condition ' conditions{cc}]);
    end
end

%% Plot group level connectivity matrices

% Make main diagonal NaN (they are really 1 but that would mess up the colormap limits
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
        for pp1 = 1:5
            for pp2 = 1:5
                count = count + 1;
                subplot(5,5,count);
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
