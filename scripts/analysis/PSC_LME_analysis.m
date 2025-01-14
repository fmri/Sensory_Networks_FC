%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to perform linear fixed effects modeling for
% percent signal change for the different task contrasts
% Tom Possidente - September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;


%% Set script variables
localizer_data = false; % localizer data instead of spacetime data
use_probabilistic_ROIs = false;
modalities_use = {'visual', 'auditory'}; % 'auditory' and/or 'visual' or 'v-a' 
domains_use = {'spatial', 'temporal'}; % temporal and/or spatial, or passive | for localizer: active or passive
use_ROItypes = false; % groups ROIs by their type (visual, auditory, MD)
include_MD = true; %  include MD ROIs or not
recruitment_variable = true; % add the recruitment variable to the table
posterior_only = false; % use only posterior ROIs
vis_aud_wm_2contrast = false; % create double contrast of visWM-visPassive - audWM-audPassive from PSCs

% reference_category = 'domain';

%% Load and format data for LME table funciton
[psc_results, perc_correct_all, hemis, modality_out, domain_out, ROIs, bad_subjs] = ...
    format_psc_data(modalities_use, domains_use, use_ROItypes, include_MD, localizer_data, use_probabilistic_ROIs);

make_categorical = false;
LME_table = create_LME_table(psc_results, perc_correct_all, hemis, modality_out, domain_out, ROIs, make_categorical, bad_subjs);

if recruitment_variable && ~localizer_data
    recruitment = ( ismember(LME_table.modality,'visual') & ismember(LME_table.domain,'temporal') ) | ( ismember(LME_table.modality,'auditory') & ismember(LME_table.domain,'spatial') );
    LME_table.recruitment = recruitment;
    LME_table.recruitment = categorical(LME_table.recruitment);
end

if posterior_only
    LME_table = LME_table(ismember(LME_table.ROItype, {'pAud','pVis'}),:);
else
    %LME_table = LME_table(~ismember(LME_table.ROItype, {'pAud','pVis'}),:);
end

if vis_aud_wm_2contrast
    LME_table = vis_aud_2contrast(LME_table);
end

%% Find smallest effect and make that the reference
% ROItypes = unique(LME_table.ROItype);
% mean_PSCs = nan(length(ROItypes),1);
% for cc = 1:length(ROItypes)
%     ROItype = ROItypes{cc};
%     in_category_mask = cell2mat(cellfun(@(x) isequal(x, ROItype), LME_table.ROItype, 'UniformOutput',false));
%     new_table = LME_table(in_category_mask,:);
%     categories = unique(new_table.(reference_category));
%     assert(length(categories)==2);
%     cat1_mask = cell2mat(cellfun(@(x) isequal(x, categories{1}), new_table.(reference_category), 'UniformOutput',false));
%     mean_PSCs(cc) = mean(new_table.PSC(cat1_mask))-mean(new_table.PSC(~cat1_mask));
% end
% [min_val, min_ind] = min(abs(mean_PSCs));
% reference = ROItypes{min_ind};
% other_inds_ordered = 1:length(ROItypes);
% other_inds_ordered = other_inds_ordered(other_inds_ordered~=min_ind);
% disp(['Reference ROI type set to ' reference ' with average PSC difference of ' num2str(min_val)]);

%LME_table.perc_correct = zscore(LME_table.perc_correct);
LME_table.subject = categorical(LME_table.subject);
LME_table.hemisphere = categorical(LME_table.hemisphere);
LME_table.domain = categorical(LME_table.domain);
LME_table.modality = categorical(LME_table.modality);
LME_table.ROItype = categorical(LME_table.ROItype);


%% Run LME model
% For passive model
% lme = fitglme(LME_table, ['PSC ~ 1 + modality * ROItype + (1 + modality + ROItype | subject) ' ...
%    ' + (1 + modality + ROItype | hemisphere)'])
%lme = fitglme(LME_table, 'PSC ~ 1 * ROItype + (1 + ROItype | subject) + (1 + ROItype | hemisphere)')

% For modality model (combined auditory and visual model)
% lme = fitglme(LME_table, ['PSC ~ 1 + modality * ROItype + (1 + modality + ROItype | subject) ' ...
%     ' + (1 + modality + ROItype | hemisphere) + (1 + modality + ROItype | perc_correct) + (1 + modality + ROItype | domain)'])
% 
% % For domain auditory or visual model
% lme = fitglme(LME_table, ['PSC ~ 1 + domain * ROItype + (1 + domain + ROItype | subject) ' ...
%    ' + (1 + domain + ROItype | hemisphere) + (1 + domain + ROItype | perc_correct) + (1 + domain + ROItype | modality)'])


% For full domain and modality model
tic
lme = fitglme(LME_table, ['PSC ~ 1 + domain * ROItype * modality + (1 + domain + ROItype + modality | subject) ' ...
   ' + (1 + domain + ROItype + modality | hemisphere) + (1 + domain + ROItype + modality | perc_correct)'])
toc

% recruitment model
% LME_table_MD = LME_table(ismember(LME_table.ROItype,'MDROI'),:);
% lme = fitglme(LME_table_MD, ['PSC ~ 1 + recruitment + (1 + recruitment | subject) ' ...
%    ' + (1 + recruitment | hemisphere) + (1 + recruitment| perc_correct)'])
% LME_table_MD = LME_table(ismember(LME_table.ROItype,{'dACC', 'preSMA', 'aINS'}),:);
% lme = fitglme(LME_table_MD, ['PSC ~ 1 + recruitment + (1 + recruitment + ROItype | subject) ' ...
%    ' + (1 + recruitment + ROItype | hemisphere) + (1 + recruitment + ROItype | perc_correct)'])

% localizer data
% lme = fitglme(LME_table, ['PSC ~ 1 + modality * ROItype + (1 + modality + ROItype | subject) ' ...
%    ' + (1 + modality + ROItype | hemisphere)'])

emm = emmeans(lme,'unbalanced');
emm.table
save('emm_lme_PSCs_norepl_modalitydomain.mat', 'emm', 'lme');
%sortrows(emm.table,'Row','descend')
plot_psc_emmeans(sortrows(emm.table,'Row','descend'));
%plot_psc_emmeans(emm.table([10,7,8,11,4,1,2,5,3,9,6],:));
%title('Spacetime | Visual Temporal WM - Auditory Temporal WM');





%%
N_cond = height(emm.table);
sigdiff_tbl = table();
for cc = 1:N_cond
    contrast = zeros(1,N_cond);
    contrast(cc) = 1;
    res_table = contrasts_wald(lme, emm, contrast);
    sigdiff_tbl = [sigdiff_tbl; {emm.table.Row{cc}, emm.table{cc,"Estimated_Marginal_Mean"}, emm.table{cc,"SE"}, res_table.pVal}];
end
sigdiff_tbl.Properties.VariableNames = {'Condition', 'EMM', 'SE', 'pVal'};

