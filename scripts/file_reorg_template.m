
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
%subjectsDir_src = [projectDir, 'data/unpacked_data_nii/'];
subjectsDir_trg = [projectDir, 'data/unpacked_data_nii/'];

%% Loop over subjs
for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};
    subjDir = [subjectsDir_trg subjCode '/bold/'];

    dir_contents = {dir(subjDir).name};
    num_runs = sum(contains(dir_contents, '00'));
    for rr = 1:num_runs
        runDir = [subjDir '00' num2str(rr) '/'];
        eventfile = [runDir '/spacetime_condition_timing.para'];
        unix(['rm ' eventfile ]);
    end
end


%ROIs
% ROI_dir = [projectDir 'data/ROIs/'];
% files = {dir(ROI_dir).name};
% files_cut = files(contains(files, 'ROIs.surf.nii'));
% 
% for ff = 1:length(files_cut)
% 
%     disp(files_cut(ff))
%     unix(['rm ' ROI_dir files_cut{ff}]);
% 
% end

