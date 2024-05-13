%%%% 
% The purpose of this script is to preprocess the T1s for the spacetime
% experiment so that they can be used in conn toolbox for surface analysis
% Created: Tom Possidente - March 2024
%%%%%

addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

unix('export SUBJECTS_DIR=/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/');

%% Set up directories and subj info

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir = [projectDir, 'data/'];

%% Start looping through subjs
for ss = 1:length(subjCodes)

    %% Set up subj varibles
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    experiment_date = subjDf_cut.([experiment_name,'Date']){subjRow};
    dirSource = [projectDir 'data/unpacked_data_nii/' subjCode, '/' ];

    assert(~contains(experiment_date,'/'), ['Subj ', subjCode, ': experiment selected has multiple dates (this should not happen)'])

    % Pick t1 date that matches experiment date if available
    t1Date = subjDf_cut.t1Date{subjRow};
    t1Run = subjDf_cut.t1Runs{subjRow};
    if contains(t1Date, '/') % if multiple dates
        t1Dates = strsplit(t1Date, '/');
        t1Runs = strsplit(t1Run, '/');
        dateMask = strcmp(t1Dates, experiment_date);
        t1Date = t1Dates{dateMask};
        t1Run = t1Runs{dateMask}; % use the one that matches experiment date

        if isempty(t1Date) 
            warning(['Subj ', subjCode, ': No t1 date matches the experiment date, using 1st t1 date available'])
            t1Date = subjDf_cut.t1Date{subjRow}{1};
            t1Run = subjDf_cut.t1Runs{subjRow}{1};
        end
    else
        if ~strcmp(t1Date, experiment_date)
            warning(['Subj ', subjCode, ': t1 date does not match the experiment date, using t1 anyway'])
        end
    end
    
    % Initialize directories 
    T1SourcePath = [dirSource 'anat/run' t1Run '/sub-' subjCode '_run' t1Run '_T1.nii'];
    
    if isfolder(['/projectnb/somerslab/scripts/jupyter/fmri/recons/' subjCode '/'])
        disp([subjCode ' already reconned, skipping.']);
    else
        continue
        %disp([subjCode ' started recon']);
        %unix(['recon-all -i ' T1SourcePath ' -subject ' subjCode ' -all']);
        %disp([subjCode ' finished recon']);
    end


end


