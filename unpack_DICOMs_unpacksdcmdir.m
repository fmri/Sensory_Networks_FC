%%%% 
% The purpose of this script is to unpack experiment-specific (trifloc space/time) DICOMs to nifti files 
% Some of this code is taken/modified from David Beeler's fmriPipeline.m function (/projectnb/somerslab/scripts/jupyter/fmri/scripts)
% Created: Tom Possidente - Feb 2024
%%%%%

addpath('/projectnb/somerslab/tom/helper_functions/')
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';
unpack_t1s = false;
unpack_func = false;
unpack_fieldmaps = true;

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
dicomsBase=[projectDir 'data/copied_DICOMs/'];

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
    dirTarget = [projectDir 'data/unpacked_data_nii/' subjCode, '/' ];
    unix(['mkdir -p ' dirTarget]); % make dir if not already there

    assert(~contains(experiment_date,'/'), ['Subj ', subjCode, ': experiment selected has multiple dates (this should not happen)'])
    
    %% Unpack t1 raw .dcm files into .nii files to anat directory
    if unpack_t1s == true
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
        t1seqName = subjDf_cut.t1SequenceName{subjRow};
        dicomsFullDir = [dicomsBase t1Date subjCode '/scans/'];
        scanSuffix = ['-', t1seqName '/resources/DICOM/files/'];
        formattedRunID = sprintf('%03d',str2double(t1Run));
        fullSrcPath = [dicomsFullDir t1Run scanSuffix];
        endTargPath = [dirTarget 'anat/sub-' subjCode '_run' t1Run '_T1.nii'];

        unix(['rm -r ' dirTarget 'anat/run' t1Run '_info'])

        if ~isfile(endTargPath) % does folder already contain .nii?
            % Actually unpack T1, then organize contents a bit
            unix(['unpacksdcmdir -src ' fullSrcPath ' -targ ' dirTarget ' -fsfast -run '...
                  t1Run ' anat nii ' 'sub-' subjCode '_run', t1Run, '_T1.nii'])
            unix(['mv ' dirTarget 'anat/' formattedRunID '/sub-' subjCode '_run', t1Run, '_T1.nii ' endTargPath]);
            unix(['mv ' dirTarget 'anat/seq.info ' dirTarget 'anat/' formattedRunID '/']);
            unix(['mv ' dirTarget 'anat/' formattedRunID ' ' dirTarget 'anat/' 'run' t1Run '_info/']);
            unix(['mv ' dirTarget  'dicomdir.sumfile ' dirTarget 'anat/' 'run' t1Run '_info/']);
            unix(['mv ' dirTarget  'unpack.log ' dirTarget 'anat/' 'run' t1Run '_info/']);
            disp([subjCode, ': T1 unpacked']);
        else
            warning(['Subj ', subjCode, ': t1 target folder already contains a .nii file. Remove this file if you want to re-unpack the t1. Skipping unpacking this t1.'])
        end
    end
    
    %% Unpack functional DICOMs to func directory as .nii
     %%% Left off here, not yet working
    % Get run numbers of func data
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);

    % Initialize directories
    func_seqName = subjDf_cut.([experiment_name, 'SequenceName']){subjRow};
    dicomsFullDir = [dicomsBase experiment_date subjCode '/scans/'];
    scanSuffix = ['-', func_seqName '/resources/DICOM/files/'];

    % Loop over runs and unpack each 
    if unpack_func
        for rr = 1:length(runs)
            run = runs(rr);
            formattedRunID = sprintf('%03d',run);
            runstr = num2str(run);
            fullSrcPath = [dicomsFullDir runstr scanSuffix];
            endTargPath = [dirTarget 'func/sub-' subjCode '_run' runstr '_' experiment_name '.nii'];
    
            if ~isfile(endTargPath) % does folder already contain .nii?
    
                %Actually unpack functional data, then organize contents a bit
                unix(['unpacksdcmdir -src ' fullSrcPath ' -targ ' dirTarget ...
                   ' -fsfast -run ' runstr ' func nii ' 'sub-' subjCode '_run' runstr '_' experiment_name '.nii'])
                unix(['mv ' dirTarget 'func/' formattedRunID '/sub-' subjCode '_run' runstr '_' experiment_name '.nii ' endTargPath]);
                unix(['mv ' dirTarget 'func/seq.info ' dirTarget 'func/' formattedRunID '/']);
                unix(['mv ' dirTarget 'func/' formattedRunID ' ' dirTarget 'func/' 'run' runstr '_info/']);
                unix(['mv ' dirTarget  'dicomdir.sumfile ' dirTarget 'func/' 'run' runstr '_info/']);
                unix(['mv ' dirTarget  'unpack.log ' dirTarget 'func/' 'run' runstr '_info/']);
                disp([subjCode ' run ' runstr ': func unpacked']);
    
            else
                warning(['Subj ' subjCode ' run ' runstr ': func target folder already contains this file. Remove this file if you want to re-unpack. Skipping.'])
            end    
        end
    end
    
    %% Unpack fieldmaps
    if unpack_fieldmaps
        % Get run numbers of func data
        FMruns = subjDf_cut.([experiment_name, 'FM']){subjRow};
        if contains(FMruns, '/') % different spacetime runs may have different fieldmaps
            FMruns = replace(FMruns, '/', ','); % still take all fieldmaps
        end
        FMruns = str2num(FMruns);        
        assert(mod(length(FMruns),2)==0, ['Subj ' subjCode 'number of fieldmaps not a multiple of 2']);

        for rr = 1:length(FMruns)
            
            run = FMruns(rr);
            formattedRunID = sprintf('%03d',run);
            runstr = num2str(run);
            if mod(rr,2) == 1
                FM_seqName = subjDf_cut.fieldmap1SequenceName{subjRow};
                fileNameSuffix = 'AP';
            else
                FM_seqName = subjDf_cut.fieldmap2SequenceName{subjRow};
                fileNameSuffix = 'PA';
            end
            scanSuffix = ['-', FM_seqName '/resources/DICOM/files/'];
            fullSrcPath = [dicomsFullDir runstr scanSuffix];

            fileName = ['/sub-' subjCode '_run' runstr '_fieldmap' fileNameSuffix '.nii'];

            if ~isfile([dirTarget 'func/' fileName]) % does folder already contain .nii?
    
                %Actually unpack fieldmap data, then organize contents a bit
                unix(['unpacksdcmdir -src ' fullSrcPath ' -targ ' dirTarget ' -fsfast -run ' runstr ' func nii ' fileName])
                unix(['mkdir -p ' dirTarget 'func/fieldmap_unpack_info/run' runstr '/'])
                unix(['rsync -a ' dirTarget 'func/' formattedRunID '// ' dirTarget '/func/fieldmap_unpack_info/run' runstr '/'])
                unix(['mv ' dirTarget 'dicomdir.sumfile ' dirTarget '/func/fieldmap_unpack_info/run' runstr '/'])
                unix(['mv ' dirTarget 'unpack.log ' dirTarget '/func/fieldmap_unpack_info/run' runstr '/'])
                unix(['mv ' dirTarget 'func/seq.info ' dirTarget '/func/fieldmap_unpack_info/run' runstr '/'])
                unix(['mv ' dirTarget 'func/fieldmap_unpack_info/run' runstr '/' fileName ' ' dirTarget 'func/'])
                unix(['rm -r ' dirTarget '/func/' formattedRunID])
            end


        end



    end

end
