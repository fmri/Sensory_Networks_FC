%%%%
% The purpose of this script is to apply the topup fieldmap conversion to
% the functional data (instead of doing this in Conn Toolbox)
% Created: Tom Possidente - May 2024
%%%%%

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

experiment_name = 'spacetime';

if strcmp(experiment_name, 'spacetime')
    fmap_prefix = '_';
    func_data_dir = 'bold';
elseif strcmp(experiment_name, 'rest')
    fmap_prefix = '_rs';
    func_data_dir = 'rest';
else
    error('Unknown experiment name');
end

%% Loop through subjs
for ss = 1:length(subjCodes)
    % Set up subj varibles
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    experiment_date = subjDf_cut.([experiment_name,'Date']){subjRow};
    if contains(experiment_date, '/') % different rs runs may have different dates
        experiment_date = strsplit(experiment_date, '/');
    else
        experiment_date = {experiment_date};
    end
    dirTarget = [projectDir 'data/unpacked_data_nii/' subjCode, '/' ];
    unix(['mkdir -p ' dirTarget]); % make dir if not already there

    % Get fieldmap run numbers
    FMruns = subjDf_cut.([experiment_name, 'FM']){subjRow};
    if contains(FMruns, '/') % different spacetime runs may have different fieldmaps
        FMruns = strsplit(FMruns, '/'); % still take all fieldmaps
    else
        FMruns = {FMruns};
    end

    for dd = 1:length(experiment_date)
        FMruns_curr = str2num(FMruns{dd});
    
        assert(mod(length(FMruns_curr),2)==0, ['Subj ' subjCode ' number of fieldmaps not a multiple of 2']);
    
        num_pairs = length(FMruns_curr)/2;
    
        % Get functional scan run numbers
        func_runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
        if contains(func_runs, '/') % runs with different fieldmaps
            func_runs = replace(func_runs, '/', ','); % still take all runs
        end
        func_runs = str2num(func_runs);
    
        % Loop through functional runs
        for ff = 1:length(func_runs)
            if strcmp(subjCode, 'NS') && strcmp(experiment_name, 'spacetime') % This subj has 1st run with 1st pair of FMs then next 3 runs with 2nd pair of FMs
                assert(num_pairs == 2, ['unexpected number of fieldmaps for subj ' subjCode ' (expected 2, got ' num2str(num_pairs)])
                if ff==1
                    topup_path = [dirTarget 'bold/sub-' subjCode 'runs' num2str(FMruns_curr(1)) num2str(FMruns_curr(2)) fmap_prefix 'fmapTopupOut'];
                else
                    topup_path = [dirTarget 'bold/sub-' subjCode 'runs' num2str(FMruns_curr(3)) num2str(FMruns_curr(4)) fmap_prefix '_fmapTopupOut'];
                end
            elseif ismember(subjCode, {'AG', 'AI'}) && strcmp(experiment_names, 'rest') % This subj has rs runs 1 and 2 with 1st pair of FMs then 3rd rs run with 2nd pair of FMs
                assert(num_pairs == 2, ['unexpected number of fieldmaps for subj ' subjCode ' (expected 2, got ' num2str(num_pairs)])
                if ff==3
                    topup_path = [dirTarget 'bold/sub-' subjCode 'runs' num2str(FMruns_curr(3)) num2str(FMruns_curr(4)) fmap_prefix 'fmapTopupOut'];
                else
                    topup_path = [dirTarget 'bold/sub-' subjCode 'runs' num2str(FMruns_curr(1)) num2str(FMruns_curr(2)) fmap_prefix 'fmapTopupOut'];
                end
            else % all other subjs should have 1 fieldmap pair
                assert(num_pairs == 1, ['unexpected number of fieldmaps for subj ' subjCode ' (expected 1, got ' num2str(num_pairs)])
                topup_path = [dirTarget 'bold/sub-' subjCode 'runs' num2str(FMruns_curr(1)) num2str(FMruns_curr(2)) fmap_prefix 'fmapTopupOut'];
            end
    
            % Apply topup to funcitonal scans
            path_func = [dirTarget func_data_dir '/00' num2str(ff) '/uf.nii'];
            path_func_out = [dirTarget func_data_dir '/00' num2str(ff) '/uf_topupApplied.nii'];
    
            if isfile(path_func_out)
                warning('There is already a functional file with topup applied in this subjects run directory, skipping...')
            else
                apply_topup_padding(path_func, path_topup_fmparams, topup_path, path_func_out, subjCode);
            end
        end
    end
    disp(['Finished subj ' subjCode]);
end


