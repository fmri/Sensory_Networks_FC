%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to apply the topup fieldmaps to
% the localizer task functional data
% Created: Tom Possidente - May 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/helper_functions/')
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/')
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
dicomsBase=[projectDir 'data/copied_DICOMs/'];
path_topup_fmparams = '/projectnb/somerslab/tom/projects/spacetime_network/data/fm_acqparams.txt';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir = [projectDir, 'data/'];

fmap_prefix = '_';
func_data_dir = 'localizer';


%% Loop through subjs
for ss = 1:length(subjCodes)
    % Set up subj varibles
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    experiment_name = 'x1WayLocalizer';
    experiment_date = subjDf_cut.([experiment_name,'Date']){subjRow};

    if isempty(experiment_date)
        experiment_name = 'x3WayLocalizer';
        experiment_date = subjDf_cut.([experiment_name,'Date']){subjRow};
    end

    if contains(experiment_date, '/')
        error(['Subj ' subjCode ' has more than one date for localizer - this should not happen'])
    else
        experiment_date = {experiment_date};
    end

    dirTarget = [projectDir 'data/unpacked_data_nii/' subjCode, '/' ];
    unix(['mkdir -p ' dirTarget]); % make dir if not already there

    % Get fieldmap run numbers
    FMruns = subjDf_cut.([experiment_name, 'FM']){subjRow};
    if contains(FMruns, '/') % different spacetime runs may have different fieldmaps
        FMruns = {replace(FMruns, '/', ',')}; % still take all fieldmaps
    else
        FMruns = {FMruns};
    end

    for dd = 1:length(experiment_date)
        FMruns_curr = str2num(FMruns{dd});

        assert(mod(length(FMruns_curr),2)==0, ['Subj ' subjCode ' number of fieldmaps not a multiple of 2']);

        % Get functional scan run numbers
        func_runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
        if contains(func_runs, '/') % runs with different fieldmaps
            func_runs = replace(func_runs, '/', ','); % still take all runs
        end
        func_runs = str2num(func_runs);

        % Loop through functional runs
        for ff = 1:length(func_runs)
            if strcmp(subjCode, 'LA') && ismember(ff,[3,4]) % functional runs 3 and 4 of subj LA use different fieldmaps
                topup_path = [dirTarget 'bold/sub-' subjCode 'runs' num2str(FMruns_curr(3)) num2str(FMruns_curr(4)) fmap_prefix 'fmapTopupOut'];
            else
                topup_path = [dirTarget 'bold/sub-' subjCode 'runs' num2str(FMruns_curr(1)) num2str(FMruns_curr(2)) fmap_prefix 'fmapTopupOut'];
            end

            % Apply topup to funcitonal scans
            path_func = [dirTarget func_data_dir '/00' num2str(ff) '/uf.nii'];
            path_func_out = [dirTarget func_data_dir '/00' num2str(ff) '/uf_topupApplied.nii'];

            if isfile(path_func_out)
                warning('There is already a functional file with topup applied in this subjects localizer directory, skipping...')
            else
                apply_topup_padding(path_func, path_topup_fmparams, topup_path, path_func_out, subjCode);
            end
        end
    end
    disp(['Finished subj ' subjCode]);
end


