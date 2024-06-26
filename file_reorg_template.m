experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir_src = [projectDir, 'data/unpacked_data_nii/'];
subjectsDir_trg = [projectDir, 'data/unpacked_data_nii_backup/'];

%% Loop over subjs and reorg dirs

for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};

    % subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    % runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
    % if contains(runs, '/') % runs with different fieldmaps
    %     runs = replace(runs, '/', ','); % still take all runs
    % end
    % runs = str2num(runs);
    % 
    % for rr = 1:length(runs)
    %     %%%%%%%%%%%%%%%%%
    %     runDir = [subjectsDir_src subjCode '/bold/00' num2str(rr) '/'];
    %     unix(['cp ' runDir '/f_events.tsv ' runDir 'sauf_topupApplied.surf_events.tsv'])
    %     %%%%%%%%%%%%%%%
    % end
end


%% ROIs
ROI_dir = [projectDir 'data/ROIs/'];
files = {dir(ROI_dir).name};
files_cut = files(contains(files, 'all_ROIs.surf.'));

for ff = 1:length(files_cut)

    disp(files_cut(ff))
    unix(['rm ' ROI_dir files_cut{ff}]);

end

