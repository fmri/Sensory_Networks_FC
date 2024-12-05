%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to preprocess the localizer data for the
%%% longDelay and Spacetime experiments using freesurfer preprocessing with
%%% topup applied fieldmaps
%%% Tom Possidente - July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
path_topup_fmparams = '/projectnb/somerslab/tom/projects/spacetime_network/data/fm_acqparams.txt';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes,{'RR', 'AH', 'SL'}));

N = length(subjCodes);

task = 'rest'; % localizer, spacetime, or rest

switch task
    case 'spacetime'
        data_dir = [projectDir 'data/unpacked_data_nii/'];
        fsd = 'bold';
        seq_name = 'spacetime';
    case 'localizer'
        data_dir = [projectDir 'data/unpacked_data_nii_fs_localizer/'];
        fsd = 'localizer';
        seq_name = 'x1WayLocalizer';
    case 'rest'
        data_dir = [projectDir 'data/unpacked_data_nii/'];
        fsd = 'rest';
        seq_name = 'rest';
end

individual_subj_space = false; % keep activation on individual subj surface (if false, use fsaverage)
use_fieldmap = true;

if use_fieldmap
    funcstem = 'fmcpr_tu';
else
    funcstem = 'fmcpr';
end

%% Loop through subjs
smooth = 3; %mm
parfor ss = 1:N
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    unix(['echo "' subjCode '" > ' data_dir subjCode '/subjectname']); % create subjectname txt file for fs

    %% Run preproc up through motion correction
    unix(['mktemplate-sess -s ' subjCode ' -d ' data_dir ' -fsd ' fsd ' -funcstem f -update']);
    unix(['register-sess -s ' subjCode ' -d ' data_dir ' -fsd ' fsd ' -delete-dat -dof 6 -bold -per-run -update -init-coreg']);
    unix(['mkbrainmask-sess -maskstem brain -fsd ' fsd ' -s ' subjCode ' -d ' data_dir ' -per-run -update']);
    unix(['mc-sess -fstem f -fmcstem fmcpr -s ' subjCode ' -d ' data_dir ' -fsd ' fsd ' -per-run -update'])

    %% Run topupapply for fieldmap distortion correction
    if use_fieldmap
        % Get fieldmap run numbers
        seq_name_curr = seq_name;
        FMruns = subjDf_cut.([seq_name_curr 'FM']){subjRow};
        if (isempty(FMruns) || strcmp(subjCode, 'MM')) && strcmp(seq_name_curr,'x1WayLocalizer')
            seq_name_curr = 'x3WayLocalizer';
            FMruns = subjDf_cut.([seq_name_curr 'FM']){subjRow};
        end
        if contains(FMruns, '/') % different spacetime runs may have different fieldmaps
            FMruns = replace(FMruns, '/', ','); % still take all fieldmaps
        end
        FMruns = str2num(FMruns);

        assert(mod(length(FMruns),2)==0, ['Subj ' subjCode ' number of fieldmaps not a multiple of 2']);

        % Get functional scan run numbers
        func_runs = subjDf_cut.([seq_name_curr, 'Runs']){subjRow};
        if contains(func_runs, '/') % runs with different fieldmaps
            func_runs = replace(func_runs, '/', ','); % still take all runs
        end
        func_runs = str2num(func_runs);

        % Loop through functional runs
        for ff = 1:length(func_runs)
            if strcmp(subjCode, 'LA') && strcmp(task, 'localizer') && ismember(ff,[3,4]) % localizer functional runs 3 and 4 of subj LA use different fieldmaps
                topup_path = [data_dir subjCode '/bold/sub-' subjCode 'runs' num2str(FMruns(3)) num2str(FMruns(4)) '_fmapTopupOut'];
            elseif strcmp(subjCode, 'NS') && strcmp(task, 'spacetime') % This subj has 1st run with 1st pair of FMs then next 3 runs with 2nd pair of FMs
                assert(length(FMruns) == 4, ['unexpected number of fieldmaps for subj ' subjCode ' (expected 4, got ' num2str(length(FMruns))])
                if ff==1
                    topup_path = [data_dir subjCode '/bold/sub-' subjCode 'runs' num2str(FMruns(1)) num2str(FMruns(2)) '_fmapTopupOut'];
                else
                    topup_path = [data_dir subjCode '/bold/sub-' subjCode 'runs' num2str(FMruns(3)) num2str(FMruns(4)) '_fmapTopupOut'];
                end
            elseif strcmp(task, 'rest')
                topup_path = [data_dir subjCode '/bold/sub-' subjCode 'runs' num2str(FMruns(1)) num2str(FMruns(2)) '_rsfmapTopupOut'];
            else
                topup_path = [data_dir subjCode '/bold/sub-' subjCode 'runs' num2str(FMruns(1)) num2str(FMruns(2)) '_fmapTopupOut'];
            end

            

            % Apply topup to funcitonal scans
            path_func = [data_dir subjCode '/' fsd '/00' num2str(ff) '/fmcpr.nii.gz'];
            path_func_out = [data_dir subjCode '/' fsd '/00' num2str(ff) '/fmcpr_tu.nii'];

            if isfile(path_func_out)
                warning('There is already a functional file with topup applied in this subjects localizer directory, skipping...')
            else
                apply_topup_padding(path_func, path_topup_fmparams, topup_path, path_func_out, subjCode);
            end
        end
    end

    %% Run the rest of preproc
    unix(['stc-sess -i ' funcstem ' -o ' funcstem '.siemens -ngroups 3 -so siemens -s ' subjCode ' -d ' data_dir ' -fsd ' fsd ' -update']);

    if individual_subj_space
        trgsubj = subjCode;
    else
        trgsubj = 'fsaverage';
    end

    unix(['rawfunc2surf-sess -i ' funcstem '.siemens -fwhm ' num2str(smooth) ' -s ' subjCode ' -d ' data_dir ' -fsd ' fsd ' -trgsubject ' trgsubj ' -stc siemens -save-unsmoothed -update -per-run']);

    disp(['Finished subj ' subjCode]);

end



