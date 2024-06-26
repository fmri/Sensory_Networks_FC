%%% The purpose of this script is to apply freesurfer motion correction
%%% to functional files
%%% Tom Possidente May 2024

%% Set up subj paths and variables
addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');

experiment_name = 'spacetime';

restingstate = true; % whether to mc resting state data as well

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
dicomsBase=[projectDir 'data/copied_DICOMs/'];
path_topup_fmparams = '/projectnb/somerslab/tom/projects/spacetime_network/data/fm_acqparams.txt';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir = [projectDir, 'data/unpacked_data_nii/'];

%% Loop through subjs
parfor ss = 1:length(subjCodes)

    %% Set up subj varibles
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);

    for rr = 1:length(runs)
        if ~isfile([subjectsDir subjCode '/bold/00' num2str(rr) '/fmcpr.nii.gz'])
            unix(['mktemplate-sess -s ' subjCode ' -d ' subjectsDir]);
            unix(['mc-sess -s ' subjCode ' -d ' subjectsDir ' -per-run']);
        else
            disp(['Motion correction file found for subj ' subjCode ' run ' num2str(rr) '... skipping']);
        end
    end

    if restingstate
        runs = subjDf_cut.('restRuns'){subjRow};
        if contains(runs, '/') % runs with different fieldmaps
            runs = replace(runs, '/', ','); % still take all runs
        end
        runs = str2num(runs);

        for rr = 1:length(runs)
            if ~isfile([subjectsDir subjCode '/rest/00' num2str(rr) '/fmcpr.nii.gz'])
                unix(['mktemplate-sess -s ' subjCode ' -d ' subjectsDir ' -fsd rest']);
                unix(['mc-sess -s ' subjCode ' -d ' subjectsDir ' -fsd rest -per-run']);
            else
                disp(['Resting state motion correction file found for subj ' subjCode ' run ' num2str(rr) '... skipping']);
            end
        end


    end
end


