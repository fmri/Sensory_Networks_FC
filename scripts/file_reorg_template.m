experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
%subjectsDir_src = [projectDir, 'data/unpacked_data_nii/'];
subjectsDir_trg = [projectDir, 'data/unpacked_data_nii_fs_localizer/'];

%% Loop over subjs
for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};

    subjDir = [subjectsDir_trg subjCode '/localizer/'];

    % add localizer contrast runlist file
    dir_contents = {dir(subjDir).name};
    num_runs = sum(contains(dir_contents, '00'));
    for rr = 1:num_runs
        files = {dir([subjDir '/00' num2str(rr) '/masks/']).name};
        files = files(contains(files, subjCode));
        for ff = 1:length(files)
            new_fname = replace(files{ff}, subjCode, 'self');
            unix(['mv ' subjDir '/00' num2str(rr) '/masks/' files{ff} ' ' subjDir '/00' num2str(rr) '/masks/' new_fname]);
        end

    end
end


    %% ROIs
    % ROI_dir = [projectDir 'data/ROIs/'];
    % files = {dir(ROI_dir).name};
    % files_cut = files(contains(files, 'all_ROIs.surf.'));
    %
    % for ff = 1:length(files_cut)
    %
    %     disp(files_cut(ff))
    %     unix(['rm ' ROI_dir files_cut{ff}]);
    %
    % end

