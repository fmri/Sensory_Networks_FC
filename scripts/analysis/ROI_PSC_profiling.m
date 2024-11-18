%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to use the PSC marginal mean/SEs for each
% task/passive conditions to visualize a functional "profile" for each ROI
% Tom Possidente - October 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Load PSC means/SEs
path_base = '/projectnb/somerslab/tom/projects/spacetime_network/data/PSCs/';

load([path_base 'posterior_passive_pscs.mat']);
post_passive_lme = lme;
post_passive_emm = emm;
post_passive_pscs = sortrows(post_passive_emm.table, "Row");

load([path_base 'passive_pscs.mat']);
passive_lme = lme;
passive_emm = emm;
passive_pscs = sortrows([passive_emm.table; post_passive_pscs], "Row");

load([path_base 'posterior_domain_modality_pscs.mat']);
post_dommod_lme = lme;
post_dommod_emm = emm;
post_dommod_pscs = sortrows(post_dommod_emm.table, "Row");

load([path_base, 'perROI_full_domain_modality_LME.mat']);
domain_modality_lme = lme;
domain_modality_emm = emmeans(domain_modality_lme,'unbalanced');
domain_modality_pscs = sortrows([domain_modality_emm.table; post_dommod_pscs], 'Row');

load([path_base, 'MD_recruitment_PSC.mat']);
recruitment_emm = emm.table;

ROIs = sort(cellstr(unique(passive_pscs.ROItype)));

%% Extract PSCs by ROI for each condition

visual_drive = passive_pscs.Estimated_Marginal_Mean(ismember(passive_pscs.modality, 'visual'));
auditory_drive = passive_pscs.Estimated_Marginal_Mean(ismember(passive_pscs.modality, 'auditory'));

visual_spatial = domain_modality_pscs.Estimated_Marginal_Mean(ismember(domain_modality_pscs.modality, 'visual') & ismember(domain_modality_pscs.domain, 'spatial'));
visual_temporal = domain_modality_pscs.Estimated_Marginal_Mean(ismember(domain_modality_pscs.modality, 'visual') & ismember(domain_modality_pscs.domain, 'temporal'));
auditory_spatial = domain_modality_pscs.Estimated_Marginal_Mean(ismember(domain_modality_pscs.modality, 'auditory') & ismember(domain_modality_pscs.domain, 'spatial'));
auditory_temporal = domain_modality_pscs.Estimated_Marginal_Mean(ismember(domain_modality_pscs.modality, 'auditory') & ismember(domain_modality_pscs.domain, 'temporal'));

%% Significance testing
% passive conditions
N_cond_passive = height(passive_emm.table);
passive_sigdiff_tbl = table();
for cc = 1:N_cond_passive
    contrast = zeros(1,N_cond_passive);
    contrast(cc) = 1;
    res_table = contrasts_wald(passive_lme, passive_emm, contrast);
    passive_sigdiff_tbl = [passive_sigdiff_tbl; {passive_emm.table.Row{cc}, passive_emm.table{cc,"Estimated_Marginal_Mean"}, passive_emm.table{cc,"SE"}, res_table.pVal}];
end
passive_sigdiff_tbl.Properties.VariableNames = {'Condition', 'EMM', 'SE', 'pVal'};
passive_sigdiff_contrast = wald_psc_emmeans(passive_lme, passive_emm);

% active conditions
N_cond = height(domain_modality_emm.table);
sigdiff_tbl = table();
for cc = 1:N_cond
    contrast = zeros(1,N_cond);
    contrast(cc) = 1;
    res_table = contrasts_wald(domain_modality_lme, domain_modality_emm, contrast);
    sigdiff_tbl = [sigdiff_tbl; {domain_modality_emm.table.Row{cc}, domain_modality_emm.table{cc,"Estimated_Marginal_Mean"}, domain_modality_emm.table{cc,"SE"}, res_table.pVal}];
end
sigdiff_tbl.Properties.VariableNames = {'Condition', 'EMM', 'SE', 'pVal'};
sigdiff_contrast = wald_psc_emmeans(domain_modality_lme, domain_modality_emm);

% passive posterior conditions
N_cond_passive_post = height(post_passive_emm.table);
postpassive_sigdiff_tbl = table();
for cc = 1:N_cond_passive_post
    contrast = zeros(1,N_cond_passive_post);
    contrast(cc) = 1;
    res_table = contrasts_wald(post_passive_lme, post_passive_emm, contrast);
    postpassive_sigdiff_tbl = [postpassive_sigdiff_tbl; {post_passive_emm.table.Row{cc}, post_passive_emm.table{cc,"Estimated_Marginal_Mean"}, post_passive_emm.table{cc,"SE"}, res_table.pVal}];
end
postpassive_sigdiff_tbl.Properties.VariableNames = {'Condition', 'EMM', 'SE', 'pVal'};
postpassive_sigdiff_contrast = wald_psc_emmeans(post_passive_lme, post_passive_emm);

% active posterior conditions
N_cond = height(post_dommod_emm.table);
sigdiff_tbl = table();
for cc = 1:N_cond
    contrast = zeros(1,N_cond);
    contrast(cc) = 1;
    res_table = contrasts_wald(post_dommod_lme, post_dommod_emm, contrast);
    sigdiff_tbl = [sigdiff_tbl; {post_dommod_emm.table.Row{cc}, post_dommod_emm.table{cc,"Estimated_Marginal_Mean"}, post_dommod_emm.table{cc,"SE"}, res_table.pVal}];
end
sigdiff_tbl.Properties.VariableNames = {'Condition', 'EMM', 'SE', 'pVal'};
sigdiff_contrast = wald_psc_emmeans(post_dommod_lme, post_dommod_emm);

%% Plot
spiderplot_data = [visual_drive, visual_spatial, auditory_spatial, auditory_drive, auditory_temporal, visual_temporal];
spiderplot_data(spiderplot_data<0) = 0;
spiderplot_data = spiderplot_data./max(spiderplot_data,[],2);
minmax = [0,1];
ax_lims = [repelem(minmax(1),6); repelem(minmax(2),6)];
colors = [repmat([0 0.4470 0.7410],4,1); repmat([0.8500 0.3250 0.0980],6,1); repmat([0.4940 0.1840 0.5560],3,1)];
colors = colors([5,6,11,7,8,12,1,2,9,3,13,4,10],:);
for rr = 1:length(ROIs)
    figure;
    s = spider_plot_class(spiderplot_data(rr,:), 'axeslimits', ax_lims, 'axesshadedlimits', {ax_lims});
    s.AxesLabels = {'visual drive', 'visual spatial', 'auditory spatial', 'auditory drive', 'auditory temporal', 'visual temporal'};
    s.FillOption = {'on'};
    s.LegendVisible = {'off'};
    s.Color = colors(rr,:);
    title(ROIs{rr});
end


