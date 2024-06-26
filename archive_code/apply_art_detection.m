%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to run SPM's ART detection on all spacetime
% subjs (post motion correction) 
% Thomas Possidente - June 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up subj paths and variables
addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
addpath('/projectnb/somerslab/tom/art_detection_spm/')

experiment_name = 'spacetime';

restingstate = true; % whether to mc resting state data as well

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
dicomsBase=[projectDir 'data/copied_DICOMs/'];
path_topup_fmparams = '/projectnb/somerslab/tom/projects/spacetime_network/data/fm_acqparams.txt';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir = [projectDir, 'data/unpacked_data_nii/'];

for ss = 1:length(subjCodes)

    %% Set up subj varibles
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);

end




