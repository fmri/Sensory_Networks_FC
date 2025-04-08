%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to loop through each subject's fs
%%% contrast results and relabel the sig.nii.gz files with the contrasts
%%% names so that they are easy to organize in freeview
%%% Tom Possidente - July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir_trg = [projectDir, 'data/unpacked_data_nii_fs_localizer/'];

%% Loop over subjs and reorg dirs
for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};

    subjDir = [subjectsDir_trg subjCode '/localizer/'];

    hem = {'lh', 'rh'};
    for hh = 1:2
        conDir = [subjDir 'localizer_contrasts_' hem{hh} '/'];
        cont_dir_contents = dir(conDir);
        subDirs = cont_dir_contents([cont_dir_contents.isdir]);
        subDirs = {subDirs(3:end).name};
        subDirs = subDirs(~ismember(subDirs, 'res'));
        for sd = 1:length(subDirs)
            unix(['cp ' conDir subDirs{sd} '/sig.nii.gz ' conDir subDirs{sd} '/' subDirs{sd} '_sig.nii.gz'])
        end
    end
    
    
end


