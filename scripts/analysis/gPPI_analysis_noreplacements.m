%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract the gPPI beta values from the
% CONN toolbox gPPI analysis results and compute group-level statistics on
% them
% Tom Possidente - October 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/'));
ccc;
%% Load missing ROIs
load('/projectnb/somerslab/tom/projects/sensory_networks_FC/data/missing_ROIs.mat', 'missing_ROIs');
load('/projectnb/somerslab/tom/projects/sensory_networks_FC/data/replacement_ROI_list.mat', 'replacement_ROIs');

%% Setup analysis parameters
plot_individual_betamaps = true;
save_out = true;

reject_subjs = {'AH', 'SL', 'RR'};
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};


task = 'vA-aA';
conditions = {'aP', 'tP', 'vP', 'aA', 'tA', 'vA', 'f'}; % Note this is a different order than the conditions in the tsv or para files, this is because Conn saves the conditions in a particular order (which is different from the order you in the condition files for some reason)

ROI_dataDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/conn_toolbox_folder/conn_localizer_task/results/firstlevel/avsm_gPPI/';
nROIs = 16*2;
missing_ROIs = []; 

skip_subjs = reject_subjs;
title_str = task;
tasks = strsplit(task,'-');
compare_conditions = cell2mat(cellfun(@(x) find(ismember(conditions, x)), tasks, 'UniformOutput', false)); % this preserves the order of the conditions
task1_perc_ind = nan;
task2_perc_ind = nan;

Nsubjs = length(subjCodes);
betas = nan(nROIs,nROIs,Nsubjs,2);

%% Load gPPI betas
load([ROI_dataDir 'resultsROI_Condition' sprintf('%03d', compare_conditions(1)) '.mat'], 'Z', 'names');
betas(:,:,:,1) = Z;
cond1_names = names;

if length(compare_conditions)==1 % if only only condition entered, assume comparing to zero
    betas(:,:,:,2) = zeros(size(Z));
else
    load([ROI_dataDir 'resultsROI_Condition' sprintf('%03d', compare_conditions(2)) '.mat'], 'Z', 'names');
    betas(:,:,:,2) = Z;
    cond2_names = names;
    assert(isequal(cond1_names, cond2_names))
end

beta_diffs = betas(:,:,:,1) - betas(:,:,:,2);
for nn = 1:length(cond1_names)
    name_preclean = strsplit(cond1_names{nn},'.');
    names{nn} = name_preclean{2}(1:end-4);
end

load('/projectnb/somerslab/tom/projects/sensory_networks_FC/data/all_replacement_ROIs.mat', 'all_replacement_ROIs');
use_replacements = false;

%% Plot mean beta matrices
mean_beta_diffs = mean(beta_diffs, 3, 'omitnan');
pVis_name = 'pVis';
desired_order = {
        'tgPCS (L)', 'FO (L)', 'CO (L)', 'cIFSG (L)', 'cmSFG (L)', 'pAud (L)', ...
        'tgPCS (R)', 'FO (R)', 'CO (R)', 'cIFSG (R)', 'cmSFG (R)', 'pAud (R)', ...
        'sPCS (L)', 'iPCS (L)', 'midIFS (L)', [pVis_name ' (L)'], 'sPCS (R)', 'iPCS (R)', 'midIFS (R)', [pVis_name ' (R)'], ...
        'sm_aINS (L)', 'sm_preSMA (L)', 'sm_dACC (L)', 'sm_sPCS (L)', 'sm_iPCS (L)', 'sm_midFSG (L)',...
        'sm_aINS (R)', 'sm_preSMA (R)', 'sm_dACC (R)', 'sm_sPCS (R)', 'sm_iPCS (R)', 'sm_midFSG (R)'};
ROI_str_mod = 5;
reject_str = {'sm_ROIs2'};
ROI_str = 'avsm_ROIs'; 
vbias_ROIs = {'sPCS', 'iPCS'};
abias_ROIs = {'tgPCS', 'cIFSG', 'cmSFG', 'CO', 'FO', 'pAud'};
sm_ROIs = {'sm_aINS', 'sm_preSMA', 'sm_dACC', 'sm_sPCS', 'sm_iPCS', 'sm_midFSG', 'midIFS'};
pVis_ROIs = {'pVis'};
pAud_ROIs = {};

for nn = 1:length(cond1_names)
    name_preclean = strsplit(cond1_names{nn},'.');
    names_clean{nn} = name_preclean{2};
end

[ROIs_match, reorder_inds] = ismember(desired_order, names_clean);
names_plot = names_clean(reorder_inds);
mean_beta_diffs = mean_beta_diffs(reorder_inds, reorder_inds);

figure;
heatmap(names_plot, names_plot, mean_beta_diffs); colormap turbo; title(task);

%% Build design matrix for LME
hemis = repelem({'lh', 'rh'},nROIs/2);

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
        ROI1 = names{rr1};
        hemi1 = hemis{rr1};
        if ismember(ROI1, abias_ROIs)
            ROI1_type = 'abias';
        elseif ismember(ROI1, vbias_ROIs)
            ROI1_type = 'vbias';
        elseif ismember(ROI1, sm_ROIs)
            ROI1_type = 'supramodal';
        elseif ismember(ROI1, pVis_ROIs)
            ROI1_type = 'pVis';
        elseif ismember(ROI1, pAud_ROIs)
            ROI1_type = 'pAud';
        else
            error('ROI type unknown');
        end
        for rr2 = 1:nROIs
            if rr1 ~= rr2 % don't include same ROI to itself
                ROI2 = names{rr2};
                hemi2 = hemis{rr2};

                if ~use_replacements
                    if ismember([subjCodes{ss} '_' ROI1 '_' hemi1], all_replacement_ROIs) | ismember([subjCodes{ss} '_' ROI2 '_' hemi2], all_replacement_ROIs)
                        continue
                    end
                end

                if ismember(ROI2, sm_ROIs)
                    ROI2_type = 'supramodal';
                elseif ismember(ROI2, abias_ROIs)
                    ROI2_type = 'abias';
                elseif ismember(ROI2, vbias_ROIs)
                    ROI2_type = 'vbias';
                elseif ismember(ROI2, pVis_ROIs)
                    ROI2_type = 'pVis';
                elseif ismember(ROI2, pAud_ROIs)
                    ROI2_type = 'pAud';
                else
                    error('ROI type unknown');
                end

                ROItype_order = sort({ROI1_type, ROI2_type});
                connection_type = [ROItype_order{1} '<->' ROItype_order{2}];
                beta_diff1 = beta_diffs(rr1,rr2,ss);
                beta_vis1 = betas(rr1,rr2,ss,1);
                beta_aud1 = betas(rr1,rr2,ss,2);

                if rr1 > rr2 % only count each occurance once
                    hemi12 = [hemi1 '_' hemi2];
                    beta_diff_mean = mean([beta_diff1, beta_diffs(rr2,rr1,ss)], 'omitnan'); % take mean of both a->b and b->a bc gPPI isn't directional and we should only use one value per connection
                    beta_vis_mean = mean([beta_vis1, betas(rr2,rr1,ss,1)], 'omitnan');
                    beta_aud_mean = mean([beta_aud1, betas(rr2,rr1,ss,2)], 'omitnan');
                    if isnan(beta_diff_mean)
                        continue;
                    end
                    data_table = [data_table; {beta_diff_mean, beta_vis_mean, beta_aud_mean, ss, hemi12, connection_type, task_pcorrect_diff}];
                end
            end
        end
    end
end

data_table.Properties.VariableNames = {'beta_diff', 'vis_beta', 'aud_beta', 'subject', 'hemispheres', 'connection_type', 'task_pcorrect_diff'};


%% Make gPPI beta bar graph

connection_types = unique(data_table.connection_type);
connection_types = connection_types(~ismember(connection_types, {'pVis<->pVis', 'pAud<->pAud'}))
n_conntypes = length(connection_types);
mean_PPI = nan(n_conntypes,2);
SE_PPI = nan(n_conntypes,2);

for cc = 1:n_conntypes
    conn_type = connection_types{cc};
    mean_PPI(cc,1) = mean(data_table.vis_beta(strcmp(data_table.connection_type, conn_type)));
    mean_PPI(cc,2) = mean(data_table.aud_beta(strcmp(data_table.connection_type, conn_type)));
    SE_PPI(cc,1) = std(data_table.vis_beta(strcmp(data_table.connection_type, conn_type)))/sqrt(sum(strcmp(data_table.connection_type, conn_type)));
    SE_PPI(cc,2) = std(data_table.aud_beta(strcmp(data_table.connection_type, conn_type)))/sqrt(sum(strcmp(data_table.connection_type, conn_type)));
end

mean_diffPPIs = abs(mean_PPI(:,1) - mean_PPI(:,2));
[~,inds] = sort(mean_diffPPIs, 'descend');
y = mean_PPI(inds,:);
err = SE_PPI(inds,:);

% Plot
figure(1); clf;
hb = bar(y); % get the bar handles
hold on;
for k = 1:size(y,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1);
end

% Set Axis properties
set(gca,'xticklabel',connection_types(inds));
ylabel('Mean PPI beta')
xlabel('Connection Type');
legend({'Visual WM', 'Auditory WM'});
title('Change in connectivity during visual and auditory working memory')
set(gca, 'FontSize', 18)

%% Change to categorical data types
data_table.subject = categorical(data_table.subject);
data_table.hemispheres = categorical(data_table.hemispheres);
data_table.connection_type = categorical(data_table.connection_type);

%% Fit LME
tic;
lme = fitglme(data_table, ['beta_diff ~ 1 + connection_type + (1 + connection_type | subject) ' ...
    ' + (1 + connection_type | hemispheres)'], ...
    'DummyVarCoding','reference');
toc


emm = emmeans(lme,'unbalanced');
emm.table
if save_out
    save(['gPPI_LME_results_localizer_' task 'avsm_noreplacements.mat'], 'lme', 'emm', 'data_table'); %%% CHANGE ME
end
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
[a,b] = sortrows(gppi_sigdiff_tbl, 'EMM', 'descend', 'ComparisonMethod', 'abs');

%% Compare Visual WM betas and Auditory WM betas 
vis_beta_conns = {'pVis<->vbias', 'pVis<->supramodal', 'supramodal<->vbias', 'vbias<->vbias'};
aud_beta_conns = {'abias<->abias', 'abias<->supramodal'};
hemisphere_possibilities = {'lh_lh', 'rh_lh', 'rh_rh'};

anova_table = table();
for ss = 1:Nsubjs
    for hh = 1:length(hemisphere_possibilities)
        vis_inds = ismember(data_table.connection_type, vis_beta_conns) & ismember(data_table.hemispheres, hemisphere_possibilities{hh}) & double(data_table.subject)==ss;
        mean_beta_vis = mean(data_table.vis_beta(vis_inds));
        anova_table = [anova_table; {mean_beta_vis, hemisphere_possibilities{hh}, ss, 1}];

        aud_inds = ismember(data_table.connection_type, aud_beta_conns) & ismember(data_table.hemispheres, hemisphere_possibilities{hh}) & double(data_table.subject)==ss;
        mean_beta_aud = mean(data_table.aud_beta(aud_inds));
        anova_table = [anova_table; {mean_beta_aud, hemisphere_possibilities{hh}, ss, 2}];
    end
end
anova_table.Properties.VariableNames = {'audvis_beta', 'hemisphere', 'subject', 'connection_type'};

[pvals,tbl,stats] = anovan(anova_table.audvis_beta, ...
    {anova_table.connection_type, anova_table.subject, anova_table.hemisphere}, ... 
    'model',1, 'random',[2,3], 'varnames',{'connection_type' 'subject', 'hemisphere'});


avg_visconns = mean(anova_table.audvis_beta(anova_table.connection_type==1));
avg_audconns = mean(anova_table.audvis_beta(anova_table.connection_type==2));
sem_visconns = std(anova_table.audvis_beta(anova_table.connection_type==1)) / sqrt(length(anova_table.audvis_beta(anova_table.connection_type==1)));
sem_audconns = std(anova_table.audvis_beta(anova_table.connection_type==2)) / sqrt(length(anova_table.audvis_beta(anova_table.connection_type==2)));

figure;
%swarmchart(anova_table.connection_type, anova_table.audvis_beta);
b1 = bar(1, avg_visconns);
hold on;
b2 = bar(2, avg_audconns);
errorbar([1;2], [avg_visconns;avg_audconns], [sem_visconns;sem_audconns], 'LineStyle','none') 
ylabel('Mean PPI Beta');
xticks([1,2])
xticklabels({'Visual Stream Connections', 'Auditory Stream Connections'});
legend([b1,b2], {'Visual WM Change in Connectivity', 'Auditory WM Change in Connectivity'});

save('gPPI_anova_plotpoints_noreplacements.mat', 'avg_visconns', 'avg_audconns', 'sem_visconns', 'sem_audconns');



%% normality checks
residuals = lme.residuals;
x = (residuals - mean(residuals))/std(residuals);
figure;
qqplot(residuals);
figure;
cdfplot(x)
hold on
x_values = linspace(min(x),max(x));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')
