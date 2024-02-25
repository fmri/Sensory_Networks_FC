%%%%
% The purpose of this script is to copy DICOMs out of BUDICOMS to a common
% directory for later unpacking. This helps ensure original DICOMs are not
% effected or accidentally edited/deleted
% Some of this code is taken/modified from David Beeler's fmriPipeline.m
% function (/projectnb/somerslab/scripts/jupyter/fmri/scripts)
%
% Note: Only the scans on the same date as the experiment_name will be
% loaded
%
% Created: Tom Possidente - Feb 2024
%%%%%

addpath('/projectnb/somerslab/tom/helper_functions/')
ccc;

% Load in subject info
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;

% Loop over subjs and copy out DICOMs
rawDicomDir = '/projectnb/somerslab/BU_DICOMs/';
projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/copied_DICOMs/';
for s=1:length(subjCodes)
    subjCode = subjCodes{s};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    subjID = [subjDf_cut.([experiment_name 'Date']){subjRow} subjCode];
    unix(['rsync -av ' rawDicomDir '/' subjID ' ' projectDir]);
end
