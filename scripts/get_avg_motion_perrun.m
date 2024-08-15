%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to plot the change in head position over time per
% subj per run
% Created: Tom Possidente - July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/')
ccc;

%%

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir_trg = [projectDir, 'data/unpacked_data_nii_fs_localizer/'];
N = length(subjCodes);

%% Loop over subjs
mean_per_run = nan(N,6);
for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};
    disp(subjCode)
    subjDir = [subjectsDir_trg subjCode '/localizer/'];

    dir_contents = {dir(subjDir).name};
    num_runs = sum(contains(dir_contents, '00'));
    for rr = 1:num_runs
        mcdata = readtable([subjDir '00' num2str(rr) '/fmcpr.mcdat'], 'FileType', 'text');
        displacement = mcdata.Var10;
        midpoint = find(displacement == 0);
        assert(length(midpoint) == 1, 'more than one timepoint with 0 motion displacement');
        displacement = displacement(~displacement==0);
        mean_per_run(ss,rr) = mean(displacement);
    end
end



