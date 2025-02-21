%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract the gPPI connectivity (correlation) values from the
% CONN toolbox gPPI analysis results and compute group-level statistics on
% them
% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Setup analysis parameters
supramodal_ROIs = false;
localizer = false;
task = 'rest';
use_replacement_ROIs = true;
plot_individual_betamaps = true;
save_out = false;


if ~localizer
    if ismember(task, {'rest', 'resting', 'rs'})
        reject_subjs = {'AH', 'RR'};
        subjCodes = {'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'AI', 'GG', 'UV', 'KQ', 'LN', 'PT', 'PL', 'NS'};
        if supramodal_ROIs
            ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_resting_state/results/firstlevel/sm_connectivity/';
        else
            ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_resting_state/results/firstlevel/connectivity_3sm/';
        end
        compare_conditions = '1';
            task1_perc_ind = nan;
        task2_perc_ind = nan;
    else
        error('non-resting state analyses not yet implemented');
    end
else
    % reject_subjs = {'AH', 'SL', 'RR', 'AI'};
    % subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS'};
    % conditions = {'aP', 'tP', 'vP', 'aA', 'tA', 'vA', 'f'}; % Note this is a different order than the conditions in the tsv or para files, this is because Conn saves the conditions in a particular order (which is different from the order you in the condition files for some reason)
    % ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_localizer_task/results/firstlevel/gPPI/';
    % title_str = task;
    % tasks = strsplit(task,'-');
    % compare_conditions = cell2mat(cellfun(@(x) find(ismember(conditions, x)), tasks, 'UniformOutput', false)); % this preserves the order of the conditions 
    % task1_perc_ind = nan;
    % task2_perc_ind = nan;
end



Nsubjs = length(subjCodes);
if supramodal_ROIs
    nROIs = 8*2;
else
    nROIs = 13*2; % 13 ROIs per hemisphere
end
corrs = nan(nROIs,nROIs,Nsubjs,2);

%% Load missing ROIs
load('/projectnb/somerslab/tom/projects/spacetime_network/data/missing_ROIs.mat', 'missing_ROIs');
load('/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/replacement_ROI_list.mat', 'replacement_ROIs');
if supramodal_ROIs
    missing_ROIs = [];
else
    missing_ROIs = [missing_ROIs'; replacement_ROIs];
end

%% Load connectivity correlation values
load([ROI_dataDir 'resultsROI_Condition00' compare_conditions '.mat'], 'Z', 'names');
corrs(:,:,:,1) = Z;
cond1_names = names;

if length(compare_conditions)==1 % if only only condition entered, assume comparing to zero
    corrs(:,:,:,2) = zeros(size(Z));
else
    load([ROI_dataDir 'resultsROI_Condition' sprintf('%03d', compare_conditions(2)) '.mat'], 'Z', 'names');
    corrs(:,:,:,2) = Z;
    cond2_names = names;
    assert(isequal(cond1_names, cond2_names))
end

corr_diffs = corrs(:,:,:,1) - corrs(:,:,:,2);
for nn = 1:length(cond1_names)
    name_preclean = strsplit(cond1_names{nn},'.');
    names{nn} = name_preclean{2}(1:end-4);
end

%% Remove missing ROIs
if ~use_replacement_ROIs
    for mm = 1:length(missing_ROIs)
        subj_ROI_hemi = strsplit(missing_ROIs{mm}, '_');
        subj_ind = find(strcmp(subj_ROI_hemi{1}, subjCodes));
        if strcmp(subj_ROI_hemi(3), 'lh')
            ROI_ind = find(ismember(names, subj_ROI_hemi(2)), 1, 'first');
        else
            ROI_ind = find(ismember(names, subj_ROI_hemi(2)), 1, 'last');
        end
        corr_diffs(ROI_ind,:,subj_ind) = NaN;
        corr_diffs(:,ROI_ind,subj_ind) = NaN;
    end
end

%% Plot mean beta matrices
mean_corr_diffs = mean(corr_diffs, 3, 'omitnan');
pVis_name = 'pVis';
if supramodal_ROIs
    desired_order = {'sm_aINS (L)', 'sm_preSMA (L)', 'sm_dACC (L)', 'sm_aINS (R)', 'sm_preSMA (R)', 'sm_dACC (R)',...
        'sm_sPCS (L)', 'sm_iPCS (L)', 'sm_midFSG (L)', 'sm_sPCS (R)', 'sm_iPCS (R)', 'sm_midFSG (R)',...
        'tgPCS (L)', 'FO (L)', 'CO (L)', 'cIFSG (L)', 'cmSFG (L)', 'tgPCS (R)', 'FO (R)', 'CO (R)', 'cIFSG (R)', 'cmSFG (R)',...
        'sPCS (L)', 'iPCS (L)', 'midIFS (L)', 'sPCS (R)', 'iPCS (R)', 'midIFS (R)',...
        'pVis (L)', 'pAud (L)', 'pVis (R)', 'pAud (R)'};
    ROI_str = 'sm_ROIs';
    reject_str = {'sm_ROIs2', 'avsm_ROIs'};
    ROI_str_mod = 5;
else
    desired_order = {'tgPCS (L)', 'FO (L)', 'CO (L)', 'cIFSG (L)', 'cmSFG (L)', 'pAud (L)', 'tgPCS (R)', 'FO (R)', 'CO (R)', ...
        'cIFSG (R)', 'cmSFG (R)', 'pAud (R)',...
        'sPCS (L)', 'iPCS (L)', 'midIFS (L)', [pVis_name ' (L)'], ...
        'sPCS (R)', 'iPCS (R)', 'midIFS (R)', [pVis_name ' (R)'], ...
        'sm_aINS (L)', 'sm_preSMA (L)', 'sm_dACC (L)', 'sm_aINS (R)', 'sm_preSMA (R)', 'sm_dACC (R)'};
    ROI_str_mod = 2;
    reject_str = {'sm_ROIs2', 'avsm_ROIs'};
    ROI_str = 'ROIs';
end

for nn = 1:length(cond1_names)
    name_preclean = strsplit(cond1_names{nn},'.');
    names_clean{nn} = name_preclean{2};
end
[ROIs_match, reorder_inds] = ismember(desired_order, names_clean);
names_plot = names_clean(reorder_inds);
mean_corr_diffs = mean_corr_diffs(reorder_inds, reorder_inds);

figure;
heatmap(names_plot, names_plot, mean_corr_diffs); colormap turbo; title(task);

%% Build design matrix for LME
vbias_ROIs = {'sPCS', 'iPCS', 'midIFS'};
abias_ROIs = {'tgPCS', 'cIFSG', 'cmSFG', 'CO', 'FO'};
mult_ROIs = {'sm_aINS', 'sm_preSMA', 'sm_dACC'};
sm_ROIs = {'sm_aINS', 'sm_preSMA', 'sm_dACC', 'sm_sPCS', 'sm_iPCS', 'sm_midFSG'};
hemis = repelem({'lh', 'rh'},13);

data_table = table();
for ss = 1:Nsubjs

    if ismember(subjCodes{ss}, reject_subjs)
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
        if supramodal_ROIs
            if ismember(ROI1, sm_ROIs)
                ROI1_type = 'supramodal';
            elseif strcmp(ROI1, 'pVis')
                ROI1_type = 'pVis';
            elseif strcmp(ROI1, 'pAud')
                ROI1_type = 'pAud';
            else
                error('ROI type unknown');
            end
        else
            if ismember(ROI1, abias_ROIs)
                ROI1_type = 'abias';
            elseif ismember(ROI1, vbias_ROIs)
                ROI1_type = 'vbias';
            elseif ismember(ROI1, mult_ROIs)
                ROI1_type = 'supramodal';
            elseif strcmp(ROI1, 'pVis')
                ROI1_type = 'pVis';
            elseif strcmp(ROI1, 'pAud')
                ROI1_type = 'pAud';
            else
                error('ROI type unknown');
            end
        end
        for rr2 = 1:nROIs
            if rr1 ~= rr2 % don't include same ROI to itself
                ROI2 = names{rr2};
                hemi2 = hemis{rr2};
                if supramodal_ROIs
                    if ismember(ROI2, sm_ROIs)
                        ROI2_type = 'supramodal';
                    elseif strcmp(ROI2, 'pVis')
                        ROI2_type = 'pVis';
                    elseif strcmp(ROI2, 'pAud')
                        ROI2_type = 'pAud';
                    else
                        error('ROI type unknown');
                    end
                else
                    if ismember(ROI2, abias_ROIs)
                        ROI2_type = 'abias';
                    elseif ismember(ROI2, vbias_ROIs)
                        ROI2_type = 'vbias';
                    elseif ismember(ROI2, mult_ROIs)
                        ROI2_type = 'supramodal';
                    elseif strcmp(ROI2, 'pVis')
                        ROI2_type = 'pVis';
                    elseif strcmp(ROI2, 'pAud')
                        ROI2_type = 'pAud';
                    else
                        error('ROI type unknown');
                    end
                end
                ROItype_order = sort({ROI1_type, ROI2_type});
                connection_type = [ROItype_order{1} '<->' ROItype_order{2}];
                beta_diff1 = corr_diffs(rr1,rr2,ss);
                if rr1 > rr2 % only count each occurance once 
                    hemi12 = [hemi1 '_' hemi2];
                    corr_diff_mean = mean([beta_diff1, corr_diffs(rr2,rr1,ss)], 'omitnan'); % take mean of both a->b and b->a bc gPPI isn't directional and we should only use one value per connection
                    if isnan(corr_diff_mean)
                        continue;
                    end
                    data_table = [data_table; {corr_diff_mean, ss, hemi12, connection_type, task_pcorrect_diff}];
                end
            end
        end
    end
end

data_table.Properties.VariableNames = {'corr_diff', 'subject', 'hemispheres', 'connection_type', 'task_pcorrect_diff'};

%% Change to categorical data types
data_table.subject = categorical(data_table.subject);
data_table.hemispheres = categorical(data_table.hemispheres);
data_table.connection_type = categorical(data_table.connection_type);

%% Fit LME
tic;
lme = fitglme(data_table, ['corr_diff ~ 1 + connection_type + (1 + connection_type | subject) ' ...
        ' + (1 + connection_type | hemispheres)'], ...
        'DummyVarCoding','reference');
toc


emm = emmeans(lme,'unbalanced');
emm.table
save(['conn_LME_results_' task '.mat'], 'lme', 'emm'); %%% CHANGE ME

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


