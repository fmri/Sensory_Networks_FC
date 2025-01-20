%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract the gPPI beta values from the
% CONN toolbox gPPI analysis results and compute group-level statistics on
% them
% Tom Possidente - October 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Setup analysis parameters
localizer = true;
plot_individual_betamaps = true;
save_out = false;

if localizer
    reject_subjs = {'AH', 'SL', 'RR', 'AI'};
    subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS'};
else
    subjCodes = {'MM'	'PP' 'MK' 'AB' 'AD'	'LA' 'AE' 'TP' 'NM'	'AF' 'AG' 'GG' 'UV'	'PQ' 'KQ' 'LN' 'RT'	'PT' 'PL' 'NS'};
    reject_subjs = {'AH', 'SL', 'RR', 'AI'};
    load('/projectnb/somerslab/tom/projects/spacetime_network/data/behavioral/behavioral/behavioral_percent_correct_data.mat', 'perc_correct_all', 'subjCodes')
    perc_correct_all = perc_correct_all(~ismember(subjCodes, reject_subjs),1:4);
    subjCodes_pcorrect = subjCodes(~ismember(subjCodes, reject_subjs));
end


task = 'vA-vP'; 
if ~localizer
    if strcmp(task, 'auditory')
        ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_auditory_allROIs/';
        compare_conditions = [6,9]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
        skip_subjs = {'LN','GG','TP'};
        title_str = 'aS - aT';
        task1_perc_ind = 3;
        task2_perc_ind = 4;
    elseif strcmp(task, 'visual')
        ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
        compare_conditions = [2,8]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
        skip_subjs = {'LA','RT'};
        title_str = 'vT - vS';
        task1_perc_ind = 1;
        task2_perc_ind = 2;
    elseif strcmp(task, 'aS-aP')
        ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
        compare_conditions = [6,3]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
        skip_subjs = {'LN','GG' ,'TP'};
        title_str = 'aS-aP';
        task1_perc_ind = nan;
        task2_perc_ind = nan;
    elseif strcmp(task, 'aT-aP')
        ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
        compare_conditions = [9,3]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
        skip_subjs = {};
        title_str = 'aT-aP';
        task1_perc_ind = nan;
        task2_perc_ind = nan;
    elseif strcmp(task, 'vT-vP')
        ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
        compare_conditions = [2,5]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
        skip_subjs = {'RT','LA'};
        title_str = 'vT-aP';
        task1_perc_ind = nan;
        task2_perc_ind = nan;
    elseif strcmp(task, 'vS-vP')
        ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
        compare_conditions = [8,5]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
        skip_subjs = {'RT','LA'};
        title_str = 'vT-aP';
        task1_perc_ind = nan;
        task2_perc_ind = nan;
    else
            error('task string not recognized');
    end
else
    conditions = {'aP', 'tP', 'vP', 'aA', 'tA', 'vA', 'f'}; % Note this is a different order than the conditions in the tsv or para files, this is because Conn saves the conditions in a particular order (which is different from the order you in the condition files for some reason)
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_localizer_task/results/firstlevel/gPPI/';
    skip_subjs = reject_subjs;
    title_str = task;
    tasks = strsplit(task,'-');
    compare_conditions = cell2mat(cellfun(@(x) find(ismember(conditions, x)), tasks, 'UniformOutput', false)); % this preserves the order of the conditions 
    task1_perc_ind = nan;
    task2_perc_ind = nan;
end



Nsubjs = length(subjCodes);
nROIs = 13*2; % 13 ROIs per hemisphere
betas = nan(nROIs,nROIs,Nsubjs,2);

%% Load gPPI betas
load([ROI_dataDir 'resultsROI_Condition00' num2str(compare_conditions(1)) '.mat'], 'Z', 'names', 'conditionnames');
if localizer
    Z = Z(:,:,1:end ~= 12); % take out index 12 which is subj AI
end
betas(:,:,:,1) = Z;
cond1_names = names;

if length(compare_conditions)==1 % if only only condition entered, assume comparing to zero
    betas(:,:,:,2) = zeros(size(Z));
else
    load([ROI_dataDir 'resultsROI_Condition00' num2str(compare_conditions(2)) '.mat'], 'Z', 'names');
    if localizer
        Z = Z(:,:,1:end ~= 12); % take out index 12 which is subj AI
    end
    betas(:,:,:,2) = Z;
    cond2_names = names;
    assert(isequal(cond1_names, cond2_names))
end

beta_diffs = betas(:,:,:,1) - betas(:,:,:,2);
names = cellfun(@(x) x(6:end-4), cond1_names, 'UniformOutput',false);

%% Plot mean beta matrices
mean_beta_diffs = mean(beta_diffs, 3);

desired_order = {'tgPCS (L)', 'FO (L)', 'CO (L)', 'cIFSG (L)', 'cmSFG (L)', 'pAud (L)', 'tgPCS (R)', 'FO (R)', 'CO (R)', ...
    'cIFSG (R)', 'cmSFG (R)', 'pAud (R)',...
    'sPCS (L)', 'iPCS (L)', 'midIFS (L)', [pVis_name ' (L)'], ...
    'sPCS (R)', 'iPCS (R)', 'midIFS (R)', [pVis_name ' (R)'], ...
    'aINS (L)', 'preSMA (L)', 'dACC (L)', 'aINS (R)', 'preSMA (R)', 'dACC (R)'};
ROI_str = 'ROIs';
names_clean = cellfun(@(f) f(length(ROI_str)+2:end), cond1_names, 'UniformOutput', false);
[ROIs_match, reorder_inds] = ismember(desired_order, names_clean);
names_plot = names_clean(reorder_inds);
mean_beta_diffs = mean_beta_diffs(reorder_inds, reorder_inds);

figure;
heatmap(names_plot, names_plot, mean_beta_diffs); colormap turbo; title(task);

%% Build design matrix for LME
vbias_ROIs = {'sPCS', 'iPCS', 'midIFS'};
abias_ROIs = {'tgPCS', 'cIFSG', 'cmSFG', 'CO', 'FO'};
mult_ROIs = {'aINS', 'preSMA', 'dACC'};
hemis = repelem({'lh', 'rh'},13);

data_table = table();
for ss = 1:Nsubjs

    if ismember(subjCodes{ss}, skip_subjs)
        continue;
    end
    if any(isnan([task1_perc_ind, task2_perc_ind])) % If comparing conditions where difficulty is irrelevant, use nan for pcorrect
        task_pcorrect_diff = nan;
    else
        task_pcorrect_diff = perc_correct_all{ss,task1_perc_ind} - perc_correct_all{ss,task2_perc_ind};
    end


    for rr1 = 1:nROIs
        ROI1 = names{rr1};
        hemi1 = hemis{rr1};
        if ismember(ROI1, abias_ROIs)
            ROI1_type = 'abias';
        elseif ismember(ROI1, vbias_ROIs)
            ROI1_type = 'vbias';
        elseif ismember(ROI1, mult_ROIs)
            ROI1_type = 'mult';
        elseif strcmp(ROI1, 'pVis')
            ROI1_type = 'pVis';
        elseif strcmp(ROI1, 'pAud')
            ROI1_type = 'pAud';
        else
            error('ROI type unknown');
        end
        for rr2 = 1:nROIs
            if rr1 ~= rr2 % don't include same ROI to itself
                ROI2 = names{rr2};
                hemi2 = hemis{rr2};
                if ismember(ROI2, abias_ROIs)
                    ROI2_type = 'abias';
                elseif ismember(ROI2, vbias_ROIs)
                    ROI2_type = 'vbias';
                elseif ismember(ROI2, mult_ROIs)
                    ROI2_type = 'mult';
                elseif strcmp(ROI2, 'pVis')
                    ROI2_type = 'pVis';
                elseif strcmp(ROI2, 'pAud')
                    ROI2_type = 'pAud';
                else
                    error('ROI type unknown');
                end
                ROItype_order = sort({ROI1_type, ROI2_type});
                connection_type = [ROItype_order{1} '<->' ROItype_order{2}];
                beta_diff1 = beta_diffs(rr1,rr2,ss);
                if rr1 > rr2 % only count each occurance once 
                    hemi12 = [hemi1 '_' hemi2];
                    beta_diff_mean = mean([beta_diff1, beta_diffs(rr2,rr1,ss)]); % take mean of both a->b and b->a bc gPPI isn't directional and we should only use one value per connection
                    data_table = [data_table; {beta_diff_mean, ss, hemi12, connection_type, task_pcorrect_diff}];
                end
            end
        end
    end
end

data_table.Properties.VariableNames = {'beta_diff', 'subject', 'hemispheres', 'connection_type', 'task_pcorrect_diff'};

%% Change to categorical data types
data_table.subject = categorical(data_table.subject);
data_table.hemispheres = categorical(data_table.hemispheres);
data_table.connection_type = categorical(data_table.connection_type);

%% Fit LME
tic;
if localizer
    lme = fitglme(data_table, ['beta_diff ~ 1 + connection_type + (1 + connection_type | subject) ' ...
        ' + (1 + connection_type | hemispheres)'], ...
        'DummyVarCoding','reference');
    % lme = fitglme(data_table, ['beta_diff ~ 1 + connection_type + (1 + connection_type | subject) ' ...
    %     ''], ...
    %     'DummyVarCoding','reference');
else
    lme = fitglme(data_table, ['beta_diff ~ 1 + connection_type + (1 + connection_type | subject) ' ...
        ' + (1 + connection_type | hemisphere1) + (1 + connection_type | hemisphere2) + (1 + connection_type | task_pcorrect_diff)'], ...
        'DummyVarCoding','reference');
end

toc


emm = emmeans(lme,'unbalanced');
emm.table
save(['gPPI_LME_results_localizer_' task '.mat'], 'lme', 'emm'); %%% CHANGE ME

plot_psc_emmeans(sortrows(emm.table,'Row','descend'));
title(title_str);

%% Sig testing
N_cond = height(emm.table);
gppi_sigdiff_tbl = table();
for cc = 1:N_cond
    contrast = zeros(1,N_cond);
    contrast(cc) = 1;
    res_table = contrasts_wald(lme, emm, contrast);
    gppi_sigdiff_tbl = [gppi_sigdiff_tbl; {emm.table.Row{cc}, emm.table{cc,"Estimated_Marginal_Mean"}, emm.table{cc,"SE"}, res_table.pVal}];
end
gppi_sigdiff_tbl.Properties.VariableNames = {'Condition', 'EMM', 'SE', 'pVal'};


