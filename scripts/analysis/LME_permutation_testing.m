%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to perform non parametric permutation
% testing on the LME models used to assess connectivity changes
% Tom Possidente - October 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%%
base_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/LME_results/';
load([base_dir 'gPPI_LME_results_localizer_vA-aA.mat'], 'emm', 'lme', 'data_table');
N = length(unique(data_table.subject));

%%
connection_inds = [2:9 11 12 14 15];
original_tstats = lme.Coefficients.tStat;
tstats_tested = original_tstats(connection_inds);
tested_names = lme.Coefficients.Name(connection_inds);

iterations = 1000;
null_dist_tstats = nan(iterations, length(tstats_tested));

tic;
parfor ii = 1:iterations
    permuted_table = table();

    % Shuffle labels for each subject
    for ss = 1:N
        subj_data = data_table(data_table.subject == categorical(ss),:);
        num_obs = height(subj_data);
        perm_inds = randperm(num_obs)';
        subj_data.connection_type = subj_data.connection_type(perm_inds);
        permuted_table = [permuted_table; subj_data];
    end
    try
        tic
        perm_lme = fitglme(permuted_table, ['beta_diff ~ 1 + connection_type + (1 + connection_type | subject) ' ...
            ' + (1 + connection_type | hemispheres)']);
        toc
        null_dist_tstats(ii,:) = perm_lme.Coefficients.tStat(connection_inds);
        parsave(['null_dist_iter' num2str(ii) '.mat'], perm_lme.Coefficients.tStat(connection_inds));
    catch e
        disp(e)
        disp('error in LME')
        parsave(['error_' num2str(ii) 'mat'],e);
    end
    disp(['Iteration ' num2str(ii) ' finished']);

end

toc
save('LME_permutation_nulldist_vA_aA.mat', "null_dist_tstats");


null_dist_tstats_max_abs = max(abs(null_dist_tstats),[],2);

pvals = (sum(tstats_tested < null_dist_tstats_max_pos)+1) / iterations;



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


