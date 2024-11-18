%%%%
% The purpose of this script is to plot percent signal change results for
% different ROIs and conditions
% Created: Tom Possidente - October 2024
%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load subject info
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.('spacetimeRuns'),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}));
N = length(subjCodes);

%% Load PSC data
subjCodes_pre = subjCodes;
load('/projectnb/somerslab/tom/projects/spacetime_network/data/PSC_results.mat', 'psc_results', 'ROIs', 'contrasts', 'subjCodes');
use_ROIs = ~ismember(ROIs,{'ppreCun'});
psc_results = psc_results(:,:,:,use_ROIs);
ROIs = ROIs(use_ROIs);
assert(all(isequal(subjCodes,subjCodes_pre))); % make sure same subjCodes in same order

N_conds = size(psc_results, 3);
N_ROIs = size(psc_results, 4);

modality_order  = {'visual', 'visual', 'auditory', 'auditory', 'visual', 'auditory'};
domain_order = {'spatial', 'temporal', 'spatial', 'temporal', 'passive', 'passive'};

ROIs_visual = {'sPCS', 'iPCS', 'midIFS'};
ROIvis_mask = ismember(ROIs, ROIs_visual)';
ROIs_auditory = {'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG'};
ROIaud_mask = ismember(ROIs, ROIs_auditory)';
ROIs_MD = {'aINS', 'dACC', 'preSMA'};
ROImd_mask = ismember(ROIs, ROIs_MD)';

psc_results = squeeze(mean(psc_results,2));

%% Plot passive condition results by ROItype
pV = squeeze(psc_results(:,5,:));
pA = squeeze(psc_results(:,6,:));
x = [1,2,3,5,6,7]; 
xlab = {'Passive Visual                            Passive Auditory'};
title_txt = 'Passive PSC';

swarmchart_PSC(x, pV, pA, [ROIvis_mask, ROIaud_mask, ROImd_mask], xlab, title_txt)

%% Plot passive condition results by ROI
N_xticks = N_ROIs*2 + N_ROIs;
xticks = 1:N_xticks;
n=3;
x1 = xticks(~ismember(xticks, xticks(n:n:end)));
title_txt = 'Passive Visual vs. Auditory PSC';
xlab = {'ROI'};
swarmchart_PSC_perROI(x1, pV, pA, ROIs, xlab, title_txt)



%% Plot modality results by ROItype
vSvT = reshape(psc_results(:,[1,2],:), N*2,N_ROIs);
aSaT = reshape(psc_results(:,[3,4],:), N*2,N_ROIs);
xlab = {'Visual Tasks                            Auditory Tasks'};
title_txt = 'Visual vs. Auditory Task PSC';
swarmchart_PSC(x, vSvT, aSaT, [ROIvis_mask, ROIaud_mask, ROImd_mask], xlab, title_txt)

%% Plot modality results by ROI
title_txt = 'Visual vs. Auditory Task PSC';
xlab = {'ROI'};
swarmchart_PSC_perROI(x1, vSvT, aSaT, ROIs, xlab, title_txt)



%% Plot domain results by ROItype
vS = squeeze(psc_results(:,1,:));
vT = squeeze(psc_results(:,2,:));
xlab = {'Visual Spatial                            Visual Temporal'};
title_txt = 'Visual Spatial vs Visual Temporal PSC';

swarmchart_PSC(x, vS, vT, [ROIvis_mask, ROIaud_mask, ROImd_mask], xlab, title_txt)

aS = squeeze(psc_results(:,3,:));
aT = squeeze(psc_results(:,4,:));
xlab = {'Auditory Spatial                            Auditory Temporal'};
title_txt = 'Auditory Spatial vs Auditory Temporal PSC';

swarmchart_PSC(x, aS, aT, [ROIvis_mask, ROIaud_mask, ROImd_mask], xlab, title_txt)

%% Plot domain results by ROI
title_txt = 'Visual Spatial vs Visual Temporal PSC';
xlab = {'ROI'};
swarmchart_PSC_perROI(x1, vS, vT, ROIs, xlab, title_txt)

title_txt = 'Auditory Spatial vs Auditory Temporal PSC';
xlab = {'ROI'};
swarmchart_PSC_perROI(x1, aS, aT, ROIs, xlab, title_txt)


%% Plot recruitment results by ROItype
vSvT = reshape(psc_results(:,[2,3],:), N*2,N_ROIs);
aSaT = reshape(psc_results(:,[1,4],:), N*2,N_ROIs);
xlab = {'Recruitment Tasks                            Non-Recruitment Tasks'};
title_txt = 'Recruitment Tasks vs. Non-Recruitment Tasks PSC';
swarmchart_PSC(x, vSvT, aSaT, [ROIvis_mask, ROIaud_mask, ROImd_mask], xlab, title_txt)

%% Plot recruitment results by ROI
title_txt = 'Recruitment Tasks vs. Non-Recruitment Tasks PSC';
xlab = {'ROI'};
swarmchart_PSC_perROI(x1, vSvT, aSaT, ROIs, xlab, title_txt)


%% 
% isROImd = cell2mat(cellfun(@(x) isequal(x, 'MDROI'), LME_table.ROItype, 'UniformOutput',false));
% isROIvis = cell2mat(cellfun(@(x) isequal(x, 'visualROI'), LME_table.ROItype, 'UniformOutput',false));
% isROIaud = cell2mat(cellfun(@(x) isequal(x, 'auditoryROI'), LME_table.ROItype, 'UniformOutput',false));
% isVis = cell2mat(cellfun(@(x) isequal(x, 'visual'), LME_table.modality, 'UniformOutput',false));
% isSpatial = cell2mat(cellfun(@(x) isequal(x, 'spatial'), LME_table.domain, 'UniformOutput',false));
% 
% LME_table.hemisphere = categorical(LME_table.hemisphere);
% LME_table.ROItype = categorical(LME_table.ROItype);
% LME_table.modality = categorical(LME_table.modality);
% LME_table.domain = categorical(LME_table.domain);
% 
% PSC_visualspatial_visROIs = mean(LME_table.PSC(isSpatial & isVis & isROIvis));
% PSC_visualspatial_audROIs = mean(LME_table.PSC(isSpatial& isVis& isROIaud));
% PSC_visualtemporal_visROIs = mean(LME_table.PSC(~isSpatial& isVis& isROIvis));
% PSC_visualtemporal_audROIs = mean(LME_table.PSC(~isSpatial& isVis& isROIaud));
% PSC_visualspatial_visROIs_se = std(LME_table.PSC(isSpatial & isVis & isROIvis))/sqrt(sum(isSpatial & isVis & isROIvis));
% PSC_visualspatial_audROIs_se = std(LME_table.PSC(isSpatial& isVis& isROIaud))/sqrt(sum(isSpatial& isVis& isROIaud));
% PSC_visualtemporal_visROIs_se = std(LME_table.PSC(~isSpatial& isVis& isROIvis))/sqrt(sum(~isSpatial & isVis & isROIvis));
% PSC_visualtemporal_audROIs_se = std(LME_table.PSC(~isSpatial& isVis& isROIaud))/sqrt(sum(~isSpatial & isVis & isROIaud));
% 
% PSC_auditoryspatial_visROIs = mean(LME_table.PSC(isSpatial& ~isVis& isROIvis));
% PSC_auditoryspatial_audROIs = mean(LME_table.PSC(isSpatial& ~isVis& isROIaud));
% PSC_auditorytemporal_visROIs = mean(LME_table.PSC(~isSpatial& ~isVis& isROIvis));
% PSC_auditorytemporal_audROIs = mean(LME_table.PSC(~isSpatial& ~isVis& isROIaud));
% PSC_auditoryspatial_visROIs_se = std(LME_table.PSC(isSpatial& ~isVis& isROIvis))/sqrt(sum(isSpatial & ~isVis & isROIvis));
% PSC_auditoryspatial_audROIs_se = std(LME_table.PSC(isSpatial& ~isVis& isROIaud))/sqrt(sum(isSpatial & ~isVis & isROIaud));
% PSC_auditorytemporal_visROIs_se = std(LME_table.PSC(~isSpatial& ~isVis& isROIvis))/sqrt(sum(~isSpatial & ~isVis & isROIvis));
% PSC_auditorytemporal_audROIs_se = std(LME_table.PSC(~isSpatial& ~isVis& isROIaud))/sqrt(sum(~isSpatial & ~isVis & isROIaud));
% 
% PSC_spatial_visROIs = mean(LME_table.PSC(isSpatial& isROIvis));
% PSC_spatial_audROIs = mean(LME_table.PSC(isSpatial& isROIaud));
% PSC_temporal_visROIs = mean(LME_table.PSC(~isSpatial& isROIvis));
% PSC_temporal_audROIs = mean(LME_table.PSC(~isSpatial& isROIaud));
% PSC_spatial_visROIs_se = std(LME_table.PSC(isSpatial& isROIvis))/sqrt(sum(isSpatial & isROIvis));
% PSC_spatial_audROIs_se = std(LME_table.PSC(isSpatial& isROIaud))/sqrt(sum(isSpatial & isROIaud));
% PSC_temporal_visROIs_se = std(LME_table.PSC(~isSpatial& isROIvis))/sqrt(sum(~isSpatial & isROIvis));
% PSC_temporal_audROIs_se = std(LME_table.PSC(~isSpatial& isROIaud))/sqrt(sum(~isSpatial & isROIaud));
% 
% PSC_spatial_mdROIs = mean(LME_table.PSC(isSpatial& isROImd));
% PSC_temporal_mdROIs = mean(LME_table.PSC(~isSpatial& isROImd));
% PSC_visualspatial_mdROIs = mean(LME_table.PSC(isSpatial& isVis& isROImd));
% PSC_visualtemporal_mdROIs = mean(LME_table.PSC(~isSpatial& isVis& isROImd));
% PSC_auditoryspatial_mdROIs = mean(LME_table.PSC(isSpatial& ~isVis& isROImd));
% PSC_auditorytemporal_mdROIs = mean(LME_table.PSC(~isSpatial& ~isVis& isROImd));
% PSC_spatial_mdROIs_se = std(LME_table.PSC(isSpatial& isROImd))/sqrt(sum(isSpatial& isROImd));
% PSC_temporal_mdROIs_se = std(LME_table.PSC(~isSpatial& isROImd))/sqrt(sum(~isSpatial& isROImd));
% PSC_visualspatial_mdROIs_se = std(LME_table.PSC(isSpatial& isVis& isROImd))/sqrt(sum(isSpatial& isVis& isROImd));
% PSC_visualtemporal_mdROIs_se = std(LME_table.PSC(~isSpatial& isVis& isROImd))/sqrt(sum(~isSpatial& isVis& isROImd));
% PSC_auditoryspatial_mdROIs_se = std(LME_table.PSC(isSpatial& ~isVis& isROImd))/sqrt(sum(isSpatial& ~isVis& isROImd));
% PSC_auditorytemporal_mdROIs_se = std(LME_table.PSC(~isSpatial& ~isVis& isROImd))/sqrt(sum(~isSpatial& ~isVis& isROImd));
% 
% %% Plot
% figure; 
% bar([1,2,3,5,6,7], [PSC_visualspatial_visROIs, PSC_visualspatial_audROIs, PSC_visualspatial_mdROIs, PSC_visualtemporal_visROIs, PSC_visualtemporal_audROIs, PSC_visualtemporal_mdROIs]);
% hold on;
% errorbar([1,2,3,5,6,7], [PSC_visualspatial_visROIs, PSC_visualspatial_audROIs, PSC_visualspatial_mdROIs, PSC_visualtemporal_visROIs, PSC_visualtemporal_audROIs, PSC_visualtemporal_mdROIs],...
%     [PSC_visualspatial_visROIs_se, PSC_visualspatial_audROIs_se, PSC_visualspatial_mdROIs_se, PSC_visualtemporal_visROIs_se, PSC_visualtemporal_audROIs_se, PSC_visualtemporal_mdROIs_se], 'LineStyle','none');
% xticklabels({'Vis ROIs', 'Aud ROIs', 'MD ROIs', 'Vis ROIs', 'Aud ROIs', 'MD ROIs'});
% xtickangle(45)
% xlabel('Spatial Task               Temporal Task');
% title('Visual Task');
% ylim([0,0.65])
% ylabel('PSC (active-passive)')
% grid;
% 
% figure;
% bar([1,2,3,5,6,7], [PSC_auditoryspatial_visROIs, PSC_auditoryspatial_audROIs, PSC_auditoryspatial_mdROIs, PSC_auditorytemporal_visROIs, PSC_auditorytemporal_audROIs, PSC_auditorytemporal_mdROIs]);
% hold on; 
% errorbar([1,2,3,5,6,7], [PSC_auditoryspatial_visROIs, PSC_auditoryspatial_audROIs, PSC_auditoryspatial_mdROIs, PSC_auditorytemporal_visROIs, PSC_auditorytemporal_audROIs, PSC_auditorytemporal_mdROIs],...
%     [PSC_auditoryspatial_visROIs_se, PSC_auditoryspatial_audROIs_se, PSC_auditoryspatial_mdROIs_se, PSC_auditorytemporal_visROIs_se, PSC_auditorytemporal_audROIs_se, PSC_auditorytemporal_mdROIs_se], 'LineStyle','none');
% xticklabels({'Vis ROIs', 'Aud ROIs', 'MD ROIs', 'Vis ROIs', 'Aud ROIs', 'MD ROIs'})
% xtickangle(45)
% xlabel('Spatial Task               Temporal Task')
% title('Auditory Task');
% ylim([0,0.65])
% ylabel('PSC (active-passive)')
% grid;
% 
% figure;
% bar([1,2,3,5,6,7], [PSC_spatial_visROIs, PSC_spatial_audROIs, PSC_spatial_mdROIs, PSC_temporal_visROIs, PSC_temporal_audROIs, PSC_temporal_mdROIs]);
% hold on;
% errorbar([1,2,3,5,6,7], [PSC_spatial_visROIs, PSC_spatial_audROIs, PSC_spatial_mdROIs, PSC_temporal_visROIs, PSC_temporal_audROIs, PSC_temporal_mdROIs],...
%     [PSC_spatial_visROIs_se, PSC_spatial_audROIs_se, PSC_spatial_mdROIs_se, PSC_temporal_visROIs_se, PSC_temporal_audROIs_se, PSC_temporal_mdROIs_se], 'LineStyle','none');
% xticklabels({'Vis ROIs', 'Aud ROIs', 'MD ROIs', 'Vis ROIs', 'Aud ROIs', 'MD ROIs'})
% xlabel('Spatial Task               Temporal Task')
% xtickangle(45)
% title('Visual and Auditory Tasks Combined');
% ylim([0,0.65])
% ylabel('PSC (active-passive)');
% grid;
% 
% 
