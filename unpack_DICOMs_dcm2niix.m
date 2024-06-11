%%%%
% The purpose of this script is to unpack experiment-specific (trifloc
% space/time) DICOMs to nifti files using dcm2niix
% Some of this code is taken/modified from David Beeler's fmriPipeline.m function (/projectnb/somerslab/scripts/jupyter/fmri/scripts)
% Created: Tom Possidente - Feb 2024
%%%%%

addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';
unpack_t1s = false;
unpack_func = false;
unpack_fieldmaps = false;
unpack_rs = true;
convert_fieldmaps = false; % convert fieldmaps from spin-echo to magnitude maps (usable by CONN)

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
        endTargPath = [dirTarget 'anat/run' t1Run '/sub-' subjCode '_run' t1Run '_T1.nii'];

        if ~isfile(endTargPath) % does folder already contain .nii?
            % Actually unpack T1
            unix(['mkdir -p ' dirTarget 'anat/run' t1Run '/']); % make dir if not already there
            cd /projectnb/somerslab/tom/helper_functions;
            unix(['dcm2niix -o ' dirTarget 'anat/run' t1Run '/' ' -f' ' sub-' subjCode '_run' ...
                t1Run '_T1 ' fullSrcPath]);
            cd /projectnb/somerslab/tom/projects/spacetime_network;
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
            endTargPath = [dirTarget 'bold/00' num2str(rr) '/f.nii'];

            if ~isfile(endTargPath) % does folder already contain .nii?
                %Actually unpack functional data
                unix(['mkdir -p ' dirTarget 'bold/00' num2str(rr) '/']); % make dir if not already there
                cd /projectnb/somerslab/tom/helper_functions;
                unix(['dcm2niix -o ' dirTarget 'bold/00' num2str(rr) '/' ' -f' ' f ' fullSrcPath]);
                cd /projectnb/somerslab/tom/projects/spacetime_network;
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

            if ~isfile([dirTarget 'bold/' fileName]) % does folder already contain .nii?
                %Actually unpack fieldmap data
                unix(['mkdir -p ' dirTarget 'bold/']); % make dir if not already there
                cd /projectnb/somerslab/tom/helper_functions;
                unix(['dcm2niix -o ' dirTarget 'bold/' ' -f' ' sub-' subjCode '_run' runstr '_fieldmap' fileNameSuffix ' ' ...
                    fullSrcPath]);
                cd /projectnb/somerslab/tom/projects/spacetime_network;
                disp([subjCode ' run ' runstr ': fieldmap unpacked']);
            else
                warning(['Subj ' subjCode ' run ' runstr ': fieldmap target folder already contains this file. Remove this file if you want to re-unpack. Skipping.'])
            end

        end
    end

    if convert_fieldmaps
        % Get run numbers of func data
        FMruns = subjDf_cut.([experiment_name, 'FM']){subjRow};
        if contains(FMruns, '/') % different spacetime runs may have different fieldmaps
            FMruns = replace(FMruns, '/', ','); % still take all fieldmaps
        end
        FMruns = str2num(FMruns);
        assert(mod(length(FMruns),2)==0, ['Subj ' subjCode 'number of fieldmaps not a multiple of 2']);

        num_pairs = length(FMruns)/2;
        run_inds = reshape(1:length(FMruns), [2,num_pairs]);

        for pp = 1:num_pairs

            fmap_runs = FMruns(run_inds(:,pp)); % get run numbers for this pair of fmaps
            fmapAP_filepath = [dirTarget 'bold/sub-' subjCode '_run' num2str(fmap_runs(1)) '_fieldmapAP.nii'];
            fmapPA_filepath = [dirTarget 'bold/sub-' subjCode '_run' num2str(fmap_runs(2)) '_fieldmapPA.nii'];
            fmapMerged_filepath = [dirTarget 'bold/sub-' subjCode 'runs' num2str(fmap_runs(1)) num2str(fmap_runs(2)) '_fmapMerged.nii.gz'];
            fmaptopup_filepath = [dirTarget 'bold/sub-' subjCode 'runs' num2str(fmap_runs(1)) num2str(fmap_runs(2)) '_fmapTopupOut'];
            %fmapMag_filepath = [dirTarget 'func/sub-' subjCode 'runs' num2str(fmap_runs(1)) num2str(fmap_runs(2)) '_fmapMag'];

            if ~isfile([fmaptopup_filepath '.nii'])
                % merge AP and PA files into one .nii
                unix(['fslmerge -t ' fmapMerged_filepath ' ' fmapAP_filepath ' ' fmapPA_filepath]);

                % Load in the merged fm and check that the dims are all even
                images = niftiread(fmapMerged_filepath);
                dims = size(images);
                odd_dims = mod(dims,2)==1;
                padded = false;
                if any(odd_dims(1:3))
                    padded = true;
                    disp(['Subj ' subjCode ': One or more dims of the in the fieldmap is odd, padding with 0s to make all dims even. Will remove padding after fieldmap transformation.']);
                    images = padarray(images, double(odd_dims), 0, 'post');
                    info = niftiinfo(fmapMerged_filepath);
                    info.Description = [info.Description ' - Used Matlab to pad dimensions'];
                    info.ImageSize = size(images);
                    info.raw.dim(2:5) = size(images);
                    split_outpath = split(fmapMerged_filepath,'.');
                    niftiwrite(images, split_outpath{1}, info);
                    unix(['rm ' fmapMerged_filepath]); % remove unpadded file
                end

                % Use topup command to convert from echo to more common fieldmaps (requires all dims of images to be even)
                unix(['topup --imain=' split_outpath{1} '.nii' ' --datain=/projectnb/somerslab/tom/projects/spacetime_network/data/fm_acqparams.txt '...
                    '--config=b02b0.cnf --iout=' fmaptopup_filepath]);

                % Use fslmaths to take time mean
                %unix(['fslmaths ' fmaptopup_filepath ' -Tmean ' fmapMag_filepath])

                if padded % if padding happened, unpad
                    images = niftiread([fmaptopup_filepath '.nii.gz']);
                    dims = size(images);
                    new_dims = dims - odd_dims; % take off 1 padded dim
                    images = images(1:new_dims(1), 1:new_dims(2), 1:new_dims(3), 1:new_dims(4));
                    info = niftiinfo([fmaptopup_filepath '.nii.gz']);
                    info.Description = [info.Description ' - Used Matlab to unpad dimensions'];
                    info.ImageSize = size(images);
                    info.raw.dim(2:5) = size(images);
                    niftiwrite(images, fmaptopup_filepath, info);
                    unix(['rm ' fmaptopup_filepath '.nii.gz']); % remove padded file

                end
            else
                warning(['Subj ' subjCode ': converted fieldmap target folder already contains this file. Remove this file if you want to re-convert. Skipping.'])
            end
        end

    end

    if unpack_rs

        runs = subjDf_cut.('restRuns'){subjRow};
        rs_date = subjDf_cut.('restDate'){subjRow};
        if isempty(rs_date)
            disp(['No resting state sessions found for subj ' subjCode '...skipping unpacking resting state DICOMs'])
            continue
        end
        if contains(runs, '/') % runs with different fieldmaps
            runs = split(runs, '/');
            dates = str2num(replace(rs_date, '/', ','));
            for ii = 1:length(runstr)
                runs{ii} = str2num(runs{ii});
            end
        else
            runs = {str2num(runs)};
            dates = str2num(rs_date);
        end

        % Initialize directories
        rs_seqName = subjDf_cut.('restSequenceName'){subjRow};
        scanSuffix = ['-', rs_seqName '/resources/DICOM/files/'];

        % Loop over runs and unpack each
        count = 0;
        for cc = 1:length(runs)
            runs_curr = runs{cc};
            dicomsFullDir = [dicomsBase num2str(dates(cc)) subjCode '/scans/'];
            for rr = 1:length(runs_curr)
                count = count + 1;
                run = runs_curr(rr);
                formattedRunID = sprintf('%03d',run);
                runstr = num2str(run);
                fullSrcPath = [dicomsFullDir runstr scanSuffix];
                endTargPath = [dirTarget 'rest/00' num2str(count) '/f.nii'];

                if ~isfile(endTargPath) % does folder already contain .nii?
                    %Actually unpack functional data
                    unix(['mkdir -p ' dirTarget 'rest/00' num2str(count) '/']); % make dir if not already there
                    cd /projectnb/somerslab/tom/helper_functions;
                    unix(['dcm2niix -o ' dirTarget 'rest/00' num2str(count) '/' ' -f' ' f ' fullSrcPath]);
                    cd /projectnb/somerslab/tom/projects/spacetime_network;
                    disp([subjCode ' run ' runstr ': resting state unpacked']);

                else
                    warning(['Subj ' subjCode ' run ' runstr ': resting state target folder already contains this file. Remove this file if you want to re-unpack. Skipping.'])
                end
            end
        end

    end


end

