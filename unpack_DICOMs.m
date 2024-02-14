%%%% 
% The purpose of this script is to unpack experiment-specific (trifloc space/time) DICOMs to nifti files 
% Some of this code is taken/modified from David Beeler's fmriPipeline.m function (/projectnb/somerslab/scripts/jupyter/fmri/scripts)
% Created: Tom Possidente - Feb 2024
%%%%%

addpath('/projectnb/somerslab/tom/helper_functions/')
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';
load_t1s = true;

projectDir = '/projectnb/somerslab/tom/projects/spatial_temp_network/';
dicomsBase='/projectnb/somerslab/scripts/jupyter/fmri/rawData/';

% force the date and t1runs columns to be read in as char instead of double (which would produce NaNs for those with multiple dates)
opts = detectImportOptions('/projectnb/somerslab/scripts/jupyter/subjectInfo.csv');
date_columns = [6,9,12,15,18,23,26,29,30,31];
vartypes = opts.VariableTypes;
vartypes(date_columns) = {'char'};
opts.VariableTypes = vartypes;

subjDf = readtable('/projectnb/somerslab/scripts/jupyter/subjectInfo.csv', opts);
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir = [projectDir, 'data/'];

%% Start looping through subjs
for ss = 1:length(subjCodes)

    %% Set up subj varibles
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    experiment_date = subjDf_cut.([experiment_name,'Date']){subjRow};

    assert(~contains(experiment_date,'/'), 'experiment selected has multiple dates (this should not happen)')
    
    %% Unpack t1 raw .dcm files into .nii files to anat directory
    if load_t1s == true
        % Pick t1 date that matches experiment date if available
        t1Date = subjDf_cut.t1Date{subjRow};
        t1Run = subjDf_cut.t1Runs{subjRow};
        if contains(t1Date, '/')
            t1Dates = strsplit(t1Date, '/');
            t1Runs = strsplit(t1Run, '/');
            dateMask = strcmp(t1Dates, experiment_date);
            t1Date = t1Dates{dateMask};
            t1Run = t1Runs{dateMask};

            if isempty(t1Date) 
                warning('No t1 date matches the experiment date, using 1st t1 date available')
                t1Date = subjDf_cut.t1Date{subjRow}{1};
                t1Run = subjDf_cut.t1Runs{subjRow}{1};
            end
        else
            if ~strcmp(t1Date, experiment_date)
                warning('t1 date does not match the experiment date, using t1 anyway')
            end
        end
                
        t1Code = subjDf_cut.t1SequenceName{subjRow};
        t1Dir_target = [projectDir 'data/' subjCode, '/' ];
        unix(['mkdir -p ' t1Dir_target]);
        dicomsFullDir=[dicomsBase t1Date subjCode '/scans/'];
        scanSuffix=['-', t1Code '/resources/DICOM/files/'];
        formattedRunID=sprintf('%03d',str2double(t1Run));

        if isempty(dir([t1Dir_target, '/anat/*.nii'])) % does folder already contain .nii?
            unix(['unpacksdcmdir -src ' dicomsFullDir t1Run scanSuffix ' -targ ' t1Dir_target ...
                ' -fsfast -run ' t1Run ' anat nii ' 'sub-' subjCode '_run', t1Run, '_T1.nii'])
            unix(['mv ' t1Dir_target 'anat/' formattedRunID '/sub-' subjCode '_run', t1Run, '_T1.nii ' t1Dir_target '/anat/'])
            unix(['mv ' t1Dir_target 'anat/seq.info ' t1Dir_target 'anat/' formattedRunID '/'])
            unix(['mv ' t1Dir_target 'anat/' formattedRunID ' ' t1Dir_target 'anat/' 'run' t1Run '_info/'])
            unix(['mv ' t1Dir_target  'dicomdir.sumfile ' t1Dir_target 'anat/' 'run' t1Run '_info/'])
            unix(['mv ' t1Dir_target  'unpack.log ' t1Dir_target 'anat/' 'run' t1Run '_info/'])

        else
            warning('t1 target folder already contains a .nii file. Remove this file if you want to re-unpack the t1. Skipping unpacking this t1.')
        end
    end
    
%     % BOLD data
%     runIDs = strsplit(subjDf_cut.([experiment_name,'Runs']){subjRow},',');
% 
%     for s=1:length(subjIDs)
%     
%         subjID=subjIDs{s}
%     
%         runIDs=str2num(runSessions{s})
%         volsDir=[projectDir '/vols_' experiment '/' subjCode '/']
%         unix(['mkdir -p ' volsDir])
%         rawDataDir=[projectDir '/copied_DICOMs/' subjID '/scans/']
%         scanSuffix=['-' sequenceName '/resources/DICOM/files/']
%         for rr=1:length(runIDs)
%             runID=runIDs(rr)
%             formattedRunID=sprintf('%03d',runID)
%             formattedRunNum=sprintf('%03d',runCounter)
%             unix(['unpacksdcmdir -src ' rawDataDir '/' num2str(runID) scanSuffix ' -targ ' volsDir ...
%                 ' -fsfast -run ' num2str(runID) ' bold nii f.nii.gz'])
%             unix(['mv ' volsDir '/bold/' formattedRunID ' ' volsDir '/bold/' formattedRunNum])
%             runCounter=runCounter+1
%         end
%     end
%     
end
