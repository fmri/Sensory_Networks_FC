%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract the connectivity (correlation) values from the
% CONN toolbox analysis results and compute group-level statistics on
% them
% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Load missing ROIs
load('/projectnb/somerslab/tom/projects/sensory_networks_FC/data/missing_ROIs.mat', 'missing_ROIs');

%% Setup analysis parameters
task = 'rest';
plot_individual_betamaps = true;
save_out = true;

if ismember(task, {'rest', 'resting', 'rs'})
    reject_subjs = {'RR','AH','PQ','RT','SL','MM'};
    subjCodes = {'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'AI', 'GG', 'UV', 'KQ', 'LN', 'PT', 'PL', 'NS'};
    compare_conditions = '1';
    task1_perc_ind = nan;
    task2_perc_ind = nan;
else
    error('non-resting state analyses not yet implemented');
end

Nsubjs = length(subjCodes);

ROI_dataDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/conn_toolbox_folder/conn_resting_state/results/firstlevel/avsm_connectivity/';
nROIs = 16*2;
missing_ROIs = [];

corrs = nan(nROIs,nROIs,Nsubjs,2);


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

%% Plot mean correlation matrices
mean_corr_diffs = mean(corr_diffs, 3, 'omitnan');
pVis_name = 'pVis';

desired_order = {
        'tgPCS (L)', 'tgPCS (R)', 'FO (L)', 'FO (R)', 'CO (L)', 'CO (R)', 'cIFSG (L)', 'cIFSG (R)', 'cmSFG (L)', 'cmSFG (R)',...
        'pAud (L)', 'pAud (R)', ...
        'sPCS (L)', 'sPCS (R)', 'iPCS (L)', 'iPCS (R)', 'midIFS (L)', 'midIFS (R)',...
        [pVis_name ' (L)'], [pVis_name ' (R)'], ...
        'sm_aINS (L)', 'sm_aINS (R)', 'sm_preSMA (L)', 'sm_preSMA (R)', 'sm_dACC (L)', 'sm_dACC (R)', 'sm_sPCS (L)', 'sm_sPCS (R)',...
        'sm_iPCS (L)', 'sm_iPCS (R)', 'sm_midFSG (L)', 'sm_midFSG (R)'};
ROI_str_mod = 5;
reject_str = {'sm_ROIs2'};
ROI_str = 'avsm_ROIs'; 
vbias_ROIs = {'sPCS', 'iPCS', 'midIFS'};
abias_ROIs = {'tgPCS', 'cIFSG', 'cmSFG', 'CO', 'FO'};
sm_ROIs = {'sm_aINS', 'sm_preSMA', 'sm_dACC', 'sm_sPCS', 'sm_iPCS', 'sm_midFSG'};
pVis_ROIs = {'pVis'};
pAud_ROIs = {'pAud'};


for nn = 1:length(cond1_names)
    name_preclean = strsplit(cond1_names{nn},'.');
    names_clean{nn} = name_preclean{2};
end
[ROIs_match, reorder_inds] = ismember(desired_order, names_clean);
names_plot = names_clean(reorder_inds);
mean_corr_diffs = mean_corr_diffs(reorder_inds, reorder_inds);

figure;
heatmap(names_plot, names_plot, mean_corr_diffs); colormap(redbluedark); title(task);
clim([-1, 1]);

%% Make Dendrogram
distance_measure = 1-abs(mean_corr_diffs);
distance_measure(1:1+nROIs:end) = 0; % Make diagonal distance zero
linkage_cluster = linkage(distance_measure, 'single');
inconsistency = inconsistent(linkage_cluster);
[clusters, cluster_txt] = linkage_output_extract(linkage_cluster, names);
figure;
[h,t,outperm] = dendrogram(linkage_cluster, nROIs, 'Labels', replace(names_plot, '_', ' '));
title('Hierarchical Clustering');
ylabel('Cluster Distance');

%% Build design matrix for LME
hemis = repelem({'lh', 'rh'},nROIs/2);
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
                corr_diff1 = corr_diffs(rr1,rr2,ss);

                if rr1 > rr2 % only count each occurance once
                    hemi12 = [hemi1 '_' hemi2];
                    corr_diff = corr_diff1;
                    if isnan(corr_diff)
                        error('corr is nan');
                    end
                    data_table = [data_table; {corr_diff, ss, hemi12, connection_type, task_pcorrect_diff}];
                end
            end
        end
    end
end

data_table.Properties.VariableNames = {'corr_diff', 'subject', 'hemispheres', 'connection_type', 'task_pcorrect_diff'};

%% Make connectivity bar graph

connection_types = unique(data_table.connection_type);
connection_types = connection_types(~ismember(connection_types, {'pVis<->pVis', 'pAud<->pAud'}));
n_conntypes = length(connection_types);
mean_conns = nan(n_conntypes,1);
SE_conns = nan(n_conntypes,1);

for cc = 1:n_conntypes
    conn_type = connection_types{cc};
    mean_conns(cc,1) = mean(data_table.corr_diff(strcmp(data_table.connection_type, conn_type)));
    SE_conns(cc,1) = std(data_table.corr_diff(strcmp(data_table.connection_type, conn_type)))/sqrt(sum(strcmp(data_table.connection_type, conn_type)));
end

[~,inds] = sort(abs(mean_conns), 'descend');
y = mean_conns(inds,:);
err = SE_conns(inds,:);
y = (exp(2.*y) - 1) ./ (exp(2.*y) + 1); % inverse fisher transform to get back to r from z

% Plot
figure(2);
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
ylim([-0.5 0.5])
ylabel('Mean Correlation')
xlabel('Connection Type');
grid on;
title('Resting State Connectivity')
set(gca, 'FontSize', 18)

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
if save_out
    save(['conn_LME_results_' task '_5networks.mat'], 'lme', 'emm'); %%% CHANGE ME
end
plot_psc_emmeans(sortrows(emm.table,'Row','descend'));

%% Sig testing
N_cond = height(emm.table);
sigdiff_tbl = table();
for cc = 1:N_cond
    contrast = zeros(1,N_cond);
    contrast(cc) = 1;
    res_table = contrasts_wald(lme, emm, contrast);
    sigdiff_tbl = [sigdiff_tbl; {emm.table.Row{cc}, emm.table{cc,"Estimated_Marginal_Mean"}, emm.table{cc,"SE"}, res_table.pVal}];
end
sigdiff_tbl.Properties.VariableNames = {'Condition', 'EMM', 'SE', 'pVal'};

%% Sig testing vbias <-> supramodal against abias <-> supramodal
res_table = contrasts_wald(lme, emm, [0 0 0 -1 0 0 0 0 0 0 0 0 0 1 0]);
res_table.pVal

figure;
b1 = bar(1, sigdiff_tbl.EMM(14));
hold on;
b2 = bar(2, sigdiff_tbl.EMM(4));
errorbar([1;2], [sigdiff_tbl.EMM(14);sigdiff_tbl.EMM(4)], [sigdiff_tbl.SE(14);sigdiff_tbl.SE(4)], 'LineStyle','none') 
ylabel('Mean Connectivity Coefficient');
xticks([1,2])
xticklabels({'frontal visual <-> supramodal', 'frontal auditory <-> supramodal'});

save('connectivity_supramodalcomp_plotpoints_5networks.mat', 'sigdiff_tbl');
