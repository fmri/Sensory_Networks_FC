%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to use the PSC marginal mean/SEs for each
% task/passive conditions to visualize a functional "profile" for each ROI
% Tom Possidente - October 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Load PSC means/SEs
path_base = '/projectnb/somerslab/tom/projects/spacetime_network/data/PSCs/';

load([path_base 'PSCs_localizer_passive_modality.mat']);
loc_passive_lme = lme;
loc_passive_emm = emm;
loc_passive_pscs = sortrows(loc_passive_emm.table, "Row");

load([path_base 'PSCs_localizer_active_modality.mat']);
loc_active_lme = lme;
loc_active_emm = emm;
loc_active_pscs = sortrows(loc_active_emm.table, "Row");

ROIs = sort(cellstr(unique(loc_passive_pscs.ROItype)));

%% Extract PSCs by ROI for each condition

visual_drive = loc_passive_pscs.Estimated_Marginal_Mean(ismember(loc_passive_pscs.modality, 'visual'));
auditory_drive = loc_passive_pscs.Estimated_Marginal_Mean(ismember(loc_passive_pscs.modality, 'auditory'));

visual_active = loc_active_pscs.Estimated_Marginal_Mean(ismember(loc_active_pscs.modality, 'visual'));
auditory_active = loc_active_pscs.Estimated_Marginal_Mean(ismember(loc_active_pscs.modality, 'auditory'));

%% Significance testing
% passive conditions
N_cond_passive = height(loc_passive_emm.table);
passive_sigdiff_tbl = table();
for cc = 1:N_cond_passive
    contrast = zeros(1,N_cond_passive);
    contrast(cc) = 1;
    res_table = contrasts_wald(loc_passive_lme, loc_passive_emm, contrast);
    passive_sigdiff_tbl = [passive_sigdiff_tbl; {loc_passive_emm.table.Row{cc}, loc_passive_emm.table{cc,"Estimated_Marginal_Mean"}, loc_passive_emm.table{cc,"SE"}, res_table.pVal}];
end
passive_sigdiff_tbl.Properties.VariableNames = {'Condition', 'EMM', 'SE', 'pVal'};
passive_sigdiff_contrast = wald_psc_emmeans(loc_passive_lme, loc_passive_emm);

% active conditions
N_cond = height(loc_active_emm.table);
sigdiff_tbl = table();
for cc = 1:N_cond
    contrast = zeros(1,N_cond);
    contrast(cc) = 1;
    res_table = contrasts_wald(loc_active_lme, loc_active_emm, contrast);
    sigdiff_tbl = [sigdiff_tbl; {loc_active_emm.table.Row{cc}, loc_active_emm.table{cc,"Estimated_Marginal_Mean"}, loc_active_emm.table{cc,"SE"}, res_table.pVal}];
end
sigdiff_tbl.Properties.VariableNames = {'Condition', 'EMM', 'SE', 'pVal'};
sigdiff_contrast = wald_psc_emmeans(loc_active_lme, loc_active_emm);


%% Plot
spiderplot_data = [visual_drive, visual_active, auditory_drive, auditory_active];
spiderplot_data(spiderplot_data<0) = 0;
spiderplot_data = spiderplot_data./max(spiderplot_data,[],2);
minmax = [0,1];
ax_lims = [repelem(minmax(1),4); repelem(minmax(2),4)];
%colors = [repmat([0 0.4470 0.7410],4,1); repmat([0.8500 0.3250 0.0980],6,1); repmat([0.4940 0.1840 0.5560],3,1)];
colors = [repmat([0 0.4470 0.7410],3,1); repmat([0.8500 0.3250 0.0980],5,1); repmat([0.4940 0.1840 0.5560],3,1)];
colors = colors([4,5,9,6,7,9,1,2,10,3,8],:);
for rr = 1:length(ROIs)
    figure;
    s = spider_plot_class(spiderplot_data(rr,:), 'axeslimits', ax_lims, 'axesshadedlimits', {ax_lims});
    s.AxesLabels = {'visual drive', 'visual active', 'auditory drive', 'auditory active'};
    s.FillOption = {'on'};
    s.LegendVisible = {'off'};
    s.Color = colors(rr,:);
    title(ROIs{rr});
end


