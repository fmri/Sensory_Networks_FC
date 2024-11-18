%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract the gPPI beta values from the
% CONN toolbox gPPI analysis results and compute group-level statistics on
% them
% Tom Possidente - October 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Setup analysis parameters
plot_individual_betamaps = true;
save_out = false;

load('/projectnb/somerslab/tom/projects/spacetime_network/data/behavioral/behavioral/behavioral_percent_correct_data.mat', 'perc_correct_all', 'subjCodes')
perc_correct_all = perc_correct_all(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}),1:4);
subjCodes_pcorrect = subjCodes(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}));

task = 'recruitment'; % auditory or visual
if strcmp(task, 'auditory')
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_auditory_allROIs/';
    compare_conditions = [6,9]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
    subj_rej = {'LN','GG','TP'};
    title_str = 'aS - aT';
    task1_perc_ind = 3;
    task2_perc_ind = 4;
elseif strcmp(task, 'visual')
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
    compare_conditions = [2,8]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
    subj_rej = {'LA','RT'};
    title_str = 'vT - vS';
    task1_perc_ind = 1;
    task2_perc_ind = 2;
elseif strcmp(task, 'aS-aP')
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
    compare_conditions = [6,3]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
    subj_rej = {'LN','GG' ,'TP'};
    title_str = 'aS-aP';
    task1_perc_ind = nan;
    task2_perc_ind = nan;
elseif strcmp(task, 'aT-aP')
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
    compare_conditions = [9,3]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
    subj_rej = {};
    title_str = 'aT-aP';
    task1_perc_ind = nan;
    task2_perc_ind = nan;
elseif strcmp(task, 'vT-vP')
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
    compare_conditions = [2,5]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
    subj_rej = {'RT','LA'};
    title_str = 'vT-aP';
    task1_perc_ind = nan;
    task2_perc_ind = nan;
elseif strcmp(task, 'vS-vP')
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
    compare_conditions = [8,5]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial, 3=auditory passive, 5=visual passive
    subj_rej = {'RT','LA'};
    title_str = 'vT-aP';
    task1_perc_ind = nan;
    task2_perc_ind = nan;
elseif strcmp(task,'recruitment')
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_allconds_allROIs/';
    compare_conditions = [6, 9; 2, 8];
    subj_rej = {{'LN','GG','TP'}; {'LA','RT'}};
    title_str = 'recruitment - non recruitment';
    task1_perc_ind = [3,1];
    task2_perc_ind = [4,2];
    perc_correct_all = [perc_correct_all; perc_correct_all];
else
    error('task string not recognized');
end

subjCodes = {'MM'	'PP' 'MK' 'AB' 'AD'	'LA' 'AE' 'TP' 'NM'	'AF' 'AG' 'GG' 'UV'	'PQ' 'KQ' 'LN' 'RT'	'PT' 'PL' 'NS'};
Nsubjs = length(subjCodes);
nROIs = 13*2; % 13 ROIs per hemisphere
recruitment_adj = length(task1_perc_ind); % this will be 2 for recruitment analysis and 1 otherwise, use as indicator

betas = nan(nROIs,nROIs,Nsubjs*recruitment_adj,2);

%% Load gPPI betas
for rr = 1:recruitment_adj
    load([ROI_dataDir 'resultsROI_Condition00' num2str(compare_conditions(rr,1)) '.mat'], 'Z', 'names');
    betas(:,:,((rr-1)*Nsubjs)+1:(Nsubjs*rr),1) = Z;
    cond1_names = names;
    load([ROI_dataDir 'resultsROI_Condition00' num2str(compare_conditions(rr,2)) '.mat'], 'Z', 'names');
    betas(:,:,((rr-1)*Nsubjs)+1:(Nsubjs*rr),2) = Z;
    cond2_names = names;
    assert(isequal(cond1_names, cond2_names))
end

beta_diffs = betas(:,:,:,1) - betas(:,:,:,2);
names = cellfun(@(x) x(6:end-4), cond1_names, 'UniformOutput',false);

%% Build design matrix for LME
vbias_ROIs = {'sPCS', 'iPCS', 'midIFS'};
abias_ROIs = {'tgPCS', 'cIFSG', 'cmSFG', 'CO', 'FO'};
mult_ROIs = {'aINS', 'preSMA', 'dACC'};
hemis = repelem({'lh', 'rh'},13);

data_table = table();
for ss = 1:Nsubjs*recruitment_adj
    if recruitment_adj==2
        indicator = (ss>Nsubjs)+1; % check visual or auditory modality
        if indicator==2
            subjID = ss-Nsubjs;
        else
            subjID = ss;
        end
        if ismember(subjCodes{subjID}, subj_rej{indicator})
            continue;
        end
        if any(isnan([task1_perc_ind(indicator), task2_perc_ind(indicator)])) % If comparing conditions where difficulty is irrelevant, use nan for pcorrect
            task_pcorrect_diff = nan;
        else
            task_pcorrect_diff = perc_correct_all{ss,task1_perc_ind(indicator)} - perc_correct_all{ss,task2_perc_ind(indicator)};
        end
    else
        if ismember(subjCodes{ss}, subj_rej)
            continue;
        end
        if any(isnan([task1_perc_ind, task2_perc_ind])) % If comparing conditions where difficulty is irrelevant, use nan for pcorrect
            task_pcorrect_diff = nan;
        else
            task_pcorrect_diff = perc_correct_all{ss,task1_perc_ind} - perc_correct_all{ss,task2_perc_ind};
        end
        subjID=ss;
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
                beta_diff = beta_diffs(rr1,rr2,ss);
                if rr1 > rr2
                    beta_diff2 = beta_diffs(rr2,rr1,ss);
                    data_table = [data_table; {mean([beta_diff, beta_diff2]), subjID, hemi1, hemi2, connection_type, task_pcorrect_diff}];
                end
            end
        end
    end
end

data_table.Properties.VariableNames = {'beta_diff', 'subject', 'hemisphere1', 'hemisphere2', 'connection_type', 'task_pcorrect_diff'};

%% Find connection type with smallest mean beta_diff (to use as reference in LME)
connection_types = unique(data_table.connection_type);
mean_betadiffs = nan(length(connection_types),1);
for cc = 1:length(connection_types)
    connection_type = connection_types{cc};
    in_category_mask = cell2mat(cellfun(@(x) isequal(x, connection_type), data_table.connection_type, 'UniformOutput',false));
    mean_betadiffs(cc) = mean(data_table.beta_diff(in_category_mask));
end
[min_val, min_ind] = min(abs(mean_betadiffs));
reference_conntype = connection_types{min_ind};
other_inds_ordered = 1:15;
other_inds_ordered = other_inds_ordered(other_inds_ordered~=min_ind);

%% Change to categorical data types
data_table.subject = categorical(data_table.subject);
data_table.hemisphere1 = categorical(data_table.hemisphere1);
data_table.hemisphere2 = categorical(data_table.hemisphere2);
data_table.connection_type = categorical(data_table.connection_type);
data_table.connection_type = reordercats(data_table.connection_type, [min_ind, other_inds_ordered]);
disp(['Reference connection type set to ' reference_conntype ' with average beta difference of ' num2str(min_val)]);

%% Fit LME
tic;
lme = fitglme(data_table, ['beta_diff ~ 1 + connection_type + (1 + connection_type | subject) ' ...
    ' + (1 + connection_type | hemisphere1) + (1 + connection_type | hemisphere2) + (1 + connection_type | task_pcorrect_diff)'], ...
    'DummyVarCoding','reference');
% lme = fitglme(data_table, ['beta_diff ~ 1 + connection_type + (1 + connection_type | subject) ' ...
%     ' + (1 + connection_type | hemisphere1) + (1 + connection_type | hemisphere2)'], ...
%     'DummyVarCoding','reference');
toc


emm = emmeans(lme,'unbalanced');
emm.table
save(['gPPI_LME_results_' task '.mat'], 'lme'); %%% CHANGE ME

plot_psc_emmeans(sortrows(emm.table,'Row','descend'));
title(title_str);
res_table = wald_psc_emmeans(lme, emm);

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


