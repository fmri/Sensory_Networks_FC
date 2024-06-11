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

addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
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

    % Copy over session with experiment data
    experimentDate = subjDf_cut.([experiment_name 'Date']){subjRow};
    subjID_exp = [experimentDate subjCode];
    unix(['rsync -av ' rawDicomDir '/' subjID_exp ' ' projectDir]);

    % Make sure we a session with the T1 data
    t1Date = subjDf_cut.t1Date{subjRow};
    t1Run = subjDf_cut.t1Runs{subjRow};
    if contains(t1Date, '/') % if multiple dates
        t1Dates = strsplit(t1Date, '/');
        t1Runs = strsplit(t1Run, '/');
        dateMask = strcmp(t1Dates, experimentDate);
        t1Date = t1Dates{dateMask};
        t1Run = t1Runs{dateMask}; % use the one that matches experiment date
    end

    if ~strcmp(t1Date, experimentDate) % if t1 date and spacetime date don't match, copy over data from t1 date as well
        unix(['rsync -av ' rawDicomDir t1Date subjCode ' ' projectDir]);
    end

    % Make sure we have all sessions with resting state
    rsDate = subjDf_cut.restDate{subjRow};
    rsRun = subjDf_cut.restRuns{subjRow};
    if contains(rsDate, '/') % if multiple dates
        rsDate = str2num(replace(rsDate, '/', ','));
        rsRun = strsplit(rsRun, '/');
    else
        rsDate = str2num(rsDate);
    end

    for dd = 1:length(rsDate)
        unix(['rsync -av ' rawDicomDir num2str(rsDate(dd)) subjCode ' ' projectDir]);
    end
end

