addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

base_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/LME_results/';
load([base_dir 'gPPI_LME_results_localizer_vA-aA.mat'], 'emm', 'lme', 'data_table');
original_tstats = lme.Coefficients.tStat;
tstats_tested = original_tstats([2:9 11 12 14 15]);
tested_names = lme.Coefficients.Name([2:9 11 12 14 15]);

base = '/projectnb/somerslab/tom/projects/spacetime_network/scripts/analysis/';
count = 0;

for ii = 1:885
    if ~isfile([base 'null_dist_iter' num2str(ii) '.mat'])
        disp(['No file for iteration ' num2str(ii)]);
        continue
    end
    count = count + 1;
    load([base 'null_dist_iter' num2str(ii) '.mat'], 'x');
    null_dist(count,:) = x;
end

save('null_dist_vA-aV.mat', 'null_dist', 'tested_names', 'tstats_tested');

%%


null_dist_tstats_max_abs = max(abs(null_dist),[],2);

for ii = 1:length(tstats_tested)
    pvals(ii) = (sum(abs(tstats_tested(ii)) < null_dist_tstats_max_abs)+1) / length(null_dist_tstats_max_abs);
end

tested_names(pvals<0.05)

