%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to
% Tom Possidente - October 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%%
base_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/';
load([base_dir 'MD_recruitment_PSC.mat'], 'emm', 'lme', 'LME_table_MD');
N = length(unique(LME_table_MD.subject));

%%
original_tstats = lme.Coefficients.tStat;
tstat_recruitment = original_tstats(2);

original_coef_names = lme.Coefficients.Name;
iterations = 5000;

null_dist_tstats = nan(iterations, 1);

tic;
for ii = 1:iterations
    permuted_table = table();

    % Shuffle recruitment labels for each subject
    for ss = 1:N
        subj_data = LME_table_MD(LME_table_MD.subject == categorical(ss),:);
        num_obs = height(subj_data);
        perm_inds = randperm(num_obs)';
        subj_data.recruitment = subj_data.recruitment(perm_inds);
        permuted_table = [permuted_table; subj_data];
    end
    try
        perm_lme = fitglme(permuted_table, ['PSC ~ 1 + recruitment + (1 + recruitment | subject) ' ...
            ' + (1 + recruitment | hemisphere) + (1 + recruitment| perc_correct)']);
        null_dist_tstats(ii,:) = perm_lme.Coefficients.tStat(2);
    catch e
        disp(e)
        disp('error in LME')
    end
end

toc

null_dist_tstats_max_pos = max(null_dist_tstats,[],2);
%null_dist_tstats_max_neg = min(null_dist_tstats,[],2);
%null_dist_tstats_max_abs = max(abs(null_dist_tstats),[],2);

pval_ROItype = (sum(tstat_recruitment < null_dist_tstats_max_pos)+1) / iterations




%%
% original_tstats = lme.Coefficients.tStat;
% tstat_ROItype = original_tstats(2);
% tstat_domain = original_tstats(3);
% tstat_ROItype_domain = original_tstats(4);
%
% original_coef_names = lme.Coefficients.Name;
% iterations = 5000;
%
% null_dist_tstats = nan(iterations, 3);
%
% tic;
% parfor ii = 1:iterations
%     permuted_table = table();
%
%     % Shuffle modality, domain, and ROI type labels for each subject
%     for ss = 1:dims(1)
%         subj_data = LME_table(LME_table.subject == ss,:);
%         num_obs = height(subj_data);
%         perm_inds = [randperm(num_obs); randperm(num_obs)]';
%         subj_data.ROItype = subj_data.ROItype(perm_inds(:,1));
%         subj_data.domain = subj_data.domain(perm_inds(:,2));
%         permuted_table = [permuted_table; subj_data];
%     end
%     try
%         perm_lme = fitlme(permuted_table, ['PSC ~ 1 + domain * ROItype + (1 + domain + ROItype | subject) ' ...
%         ' + (1 + domain + ROItype | hemisphere) + (1 + domain + ROItype | perc_correct)']);
%         null_dist_tstats(ii,:) = perm_lme.Coefficients.tStat(2:4);
%     catch e
%         disp(e)
%         disp('error in LME')
%     end
% end
%
% toc
%
% null_dist_tstats_max_pos = max(null_dist_tstats,[],2);
% null_dist_tstats_max_neg = min(null_dist_tstats,[],2);
% null_dist_tstats_max_abs = max(abs(null_dist_tstats),[],2);
%
% pval_ROItype = (sum(tstat_ROItype < null_dist_tstats_max_abs)+1) / iterations
% pval_domain = (sum(abs(tstat_domain) < null_dist_tstats_max_abs)+1) / iterations
% pval_ROItype_domain = (sum(tstat_ROItype_domain > null_dist_tstats_max_neg)+1) / iterations


