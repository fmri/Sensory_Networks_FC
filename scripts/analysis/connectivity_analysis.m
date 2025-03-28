%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract the ROI timeseries data for each
% subject and use correlation analysis to calculate functional connectivity
% as well as hierarchical clustering to do network analysis
% Tom Possidente - June 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Setup analysis parameters
use_data = 'rest'; % 'rest', 'localizer', or 'spacetime'
use_replacements = true; % whether to use or skip 15% replacement ROIs
ROI_set = 3; % 1 or 2 o 3
hierarchical_clustering = true;
plot_individual_connmats = false;
save_out = false;
bootstrap_hca = true;
bootstrap_iters = 10000;

switch use_data
    case 'rest'
        ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_resting_state/results/preprocessing/';
        subjCodes = {'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'AI', 'GG', 'UV', 'KQ', 'LN', 'PT', 'PL', 'NS'};
        conditions = {'rest'};
        reject_conditions = {{}};
    case 'localizer'
        ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_localizer_task/results/preprocessing/';
        subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS'};
        conditions = {'aP', 'tP', 'vP', 'aA', 'tA', 'vA', 'f'}; %Note this is a different order than the conditions in the tsv or para files, this is because Conn saves the conditions in a particular order (which is different from the order you in the condition files for some reason)
        ROI_str = 'ROIs'; % "ROI_mod" for pVis without DO, "ROI" for pVis with DO
        reject_conditions = {{},{},{},{},{},{},{}};
    case 'spacetime'
        ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/preprocessing/';
        subjCodes = {'MM'	'PP' 'MK' 'AB' 'AD'	'LA' 'AE' 'TP' 'NM'	'AF' 'AG' 'GG' 'UV'	'PQ' 'KQ' 'LN' 'RT'	'PT' 'PL' 'NS'};
        conditions = {'f' , 'aP', 'tP', 'vP', 'aS', 'tS', 'vS', 'aT', 'tT', 'vT'}; % I don't think this is the correct order of conditions - conn saves the order weird so double check
        ROI_str = 'ROIs';
        reject_conditions = {{}, {}, {}, {}, {'LN', 'GG', 'TP'}, {}, {}, {'RT', 'LA'}, {}, {}}; % These subjs had below 55% accuracy on these tasks and should be rejected from analysis in those conditions
end

%% Setup ROIs
pVis_name = 'pVis'; % pVis_mod for pVis wihout DO, pVis for pVis with DO

if ROI_set == 1
    desired_order = {'tgPCS (L)', 'FO (L)', 'CO (L)', 'cIFSG (L)', 'cmSFG (L)', 'pAud (L)', 'tgPCS (R)', 'FO (R)', 'CO (R)', ...
        'cIFSG (R)', 'cmSFG (R)', 'pAud (R)',...
        'sPCS (L)', 'iPCS (L)', 'midIFS (L)', [pVis_name ' (L)'], ...
        'sPCS (R)', 'iPCS (R)', 'midIFS (R)', [pVis_name ' (R)'], ...
        'sm_aINS (L)', 'sm_preSMA (L)', 'sm_dACC (L)', 'sm_aINS (R)', 'sm_preSMA (R)', 'sm_dACC (R)'};
    ROI_str_mod = 2;
    reject_str = {'sm_ROIs2', 'avsm_ROIs'};
    ROI_str = 'ROIs'; 

elseif ROI_set == 2
    desired_order = {'pAud (L)', 'pAud (R)',...
        [pVis_name ' (L)'], [pVis_name ' (R)'], ...
        'sm_aINS (L)', 'sm_preSMA (L)', 'sm_dACC (L)', 'sm_sPCS (L)', 'sm_iPCS (L)', 'sm_midFSG (L)',...
        'sm_aINS (R)', 'sm_preSMA (R)', 'sm_dACC (R)', 'sm_sPCS (R)', 'sm_iPCS (R)', 'sm_midFSG (R)'};
    ROI_str_mod = 5;
    reject_str = {'sm_ROIs2', 'avsm_ROIs'};
    ROI_str = 'ROIs'; %
elseif ROI_set == 3
    desired_order = {
        'tgPCS (L)', 'FO (L)', 'CO (L)', 'cIFSG (L)', 'cmSFG (L)', 'pAud (L)', ...
        'tgPCS (R)', 'FO (R)', 'CO (R)', 'cIFSG (R)', 'cmSFG (R)', 'pAud (R)', ...
        'sPCS (L)', 'iPCS (L)', 'midIFS (L)', [pVis_name ' (L)'], 'sPCS (R)', 'iPCS (R)', 'midIFS (R)', [pVis_name ' (R)'], ...
        'sm_aINS (L)', 'sm_preSMA (L)', 'sm_dACC (L)', 'sm_sPCS (L)', 'sm_iPCS (L)', 'sm_midFSG (L)',...
        'sm_aINS (R)', 'sm_preSMA (R)', 'sm_dACC (R)', 'sm_sPCS (R)', 'sm_iPCS (R)', 'sm_midFSG (R)'};
    ROI_str_mod = 5;
    reject_str = {'sm_ROIs2'};
    ROI_str = 'avsm_ROIs'; 
else
    error('ROI set beyond 1 or 2 have not been set up');
end

N = length(subjCodes);
Ncond = length(conditions);
ROI_data = cell(N,Ncond);
first_ROI_check = [ROI_str '.CO (L)'];


%% Get missing ROI data
load('/projectnb/somerslab/tom/projects/spacetime_network/data/missing_ROIs.mat', 'missing_ROIs');
load('/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/replacement_ROI_list.mat', 'replacement_ROIs');
missing_ROIs = [missing_ROIs'; replacement_ROIs];

%% Get ROI data

files = {dir(ROI_dataDir).name};

for ss = 1:N
    subjnum = sprintf( '%03d', ss ) ;
    subjfiles = files(contains(files, ['ROI_Subject' subjnum '_Condition']) & ~contains(files, '000'));

    condcount = 0;
    for ff = 1:length(subjfiles)

        load([ROI_dataDir subjfiles{ff}], 'names', 'data', 'conditionname', 'conditionweights');
        ROI_names_mask = cellfun(@(x) contains(x,ROI_str) & ~contains(x, reject_str), names);
        names = names(ROI_names_mask);
        data = data(ROI_names_mask);

        switch use_data
            case 'rest'
                assert(strcmp(conditionname,'rest'), ['Condition name is not "rest" for subj ' num2str(ss)])
                condcount = condcount + 1;
                cond_ind = 1;
            case {'localizer', 'spacetime'}
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
        for nn = 1:length(names)
            name_preclean = strsplit(names{nn},'.');
            names_clean{nn} = name_preclean{2};
        end
        

        % Reorder ROIs
        [ROIs_match, reorder_inds] = ismember(desired_order, names_clean);
        names = names_clean(reorder_inds);
        data = data(reorder_inds);
        ROI_data{ss,cond_ind} = data;

    end

    switch use_data
        case 'rest'
            assert(condcount==1, ['Unexpected number of conditions for subj ' num2str(ss)]);
        case 'localizer'
            assert(condcount==7, ['Unexpected number of conditions for subj ' num2str(ss)]);
        case 'spacetime'
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
                if (ismember(ROI2_name, missing_ROIs) || ismember(ROI1_name, missing_ROIs)) && ~use_replacements % if subj is missing ROI, connmat is nan
                    connmats(rr1,rr2,ss,cc) = NaN;
                    pvals(rr1,rr2,ss,cc) = NaN;
                    disp(['Subj ' num2str(ss) ' (' subjID ') missing ROIs: ' missing_ROIs{ismember(missing_ROIs,ROI2_name)} missing_ROIs{ismember(missing_ROIs,ROI1_name)}]);
                elseif rr2==rr1 || isempty(ROI1) || isempty(ROI2) % skip if assessing self-connectivity, or if there is no data for that condition in that subj and ROI
                    continue
                else
                    [connmats(rr1,rr2,ss,cc), pvals(rr1,rr2,ss,cc)] = corr(ROI1, ROI2);
                end
            end
        end
    end
end

%% Connectivity group
% convert to Z scores here before averaging, then convert back to corrs (r)
connmats_z = atanh(connmats);
connmat_group_z = squeeze(mean(connmats_z, 3, 'omitnan'));
connmat_group = (exp(2.*connmat_group_z) - 1) ./ (exp(2.*connmat_group_z) + 1); % inverse fisher transform

%% Calculate hierarchical clustering
if hierarchical_clustering
    for cc = 1:Ncond
        distance_measure = 1-abs(connmat_group(:,:,cc));
        distance_measure(1:1+N_ROIs:end) = 0; % Make diagonal distance zero
        linkage_cluster = linkage(distance_measure, 'ward');
        inconsistency = inconsistent(linkage_cluster);
        [clusters, cluster_txt] = linkage_output_extract(linkage_cluster, names);

        cluster_stats = HCA_optimal_cluster(clusters, distance_measure, names);
        figure; plot(cluster_stats.Num_clusters, cluster_stats.within_cluster_sum);
        grid on; xlabel('Number of Clusters'); ylabel('Within-Cluster Sum of Distance');
        if bootstrap_hca
            tic;
            cluster_prob = HCA_bootstrap(connmats(:,:,:,cc), bootstrap_iters, clusters, cluster_txt, names);
            toc;
        end
        figure;
        [h,t,outperm] = dendrogram(linkage_cluster, N_ROIs, 'Labels', replace(names, '_', ' '));
        title(['Hierarchical Clustering Condition ' conditions{cc}]);
        ylabel('Cluster Distance');
    end
end

%% Plot group level connectivity matrices

% Make main diagonal NaN (they are really 1 but that would mess up the colormap limits)
connmat_group_nandiag = NaN(size(connmat_group));
for mm = 1:size(connmat_group,3)
    connmat_group_nandiag(:,:,mm) = connmat_group(:,:,mm) - diag(diag(connmat_group(:,:,mm))) + diag(NaN(size(connmat_group,1),1));
end

for cc = 1:Ncond
    figure;
    heatmap(replace(names, '_', ' '), replace(names, '_', ' '), connmat_group_nandiag(:,:,cc), 'Colormap', turbo, 'ColorLimits', [min(connmat_group_nandiag(),[], 'all'),max(connmat_group_nandiag,[],'all')], 'FontSize',16)
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
switch use_data
    case 'localizer'
        %ztrans_connmat_group = atanh(connmat_group_nandiag);% fisher z-transoform the rs
        diffs1 = [6,4, 6, 3,1];
        diffs2 = [3,1, 4, 7,7];   % {'aP', 'tP', 'vP', 'aA', 'tA', 'vA', 'f'};
        for cc = 1:length(diffs1)
            figure;
            diff = connmat_group_nandiag(:,:,diffs1(cc))-connmat_group_nandiag(:,:,diffs2(cc));
            %diff = (exp(2.*diff_z) - 1) ./ (exp(2.*diff_z) + 1); % inverse fisher transform
            heatmap(names, names,  diff, 'Colormap', turbo, 'ColorLimits', [min(diff(),[], 'all'),max(diff,[],'all')], 'FontSize',16)
            title([conditions{diffs1(cc)} '-' conditions{diffs2(cc)}])
        end
    case 'spacetime'
        ztrans_connmat_group = atanh(connmat_group_nandiag);% fisher z-transoform the rs
        diffs1 = [2,4, 5,8,  7,10, 5,7];
        diffs2 = [1,1, 2,2 , 4,4,  8,10];
        for cc = 1:length(diffs1)
            figure;
            diff_z = ztrans_connmat_group(:,:,diffs1(cc))-ztrans_connmat_group(:,:,diffs2(cc));
            diff = (exp(2.*diff_z) - 1) ./ (exp(2.*diff_z) + 1); % inverse fisher transform
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
