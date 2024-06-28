%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to run topup on the resting state fieldmaps
% and then apply them to the functional resting state data.
% Created: Tom Possidente - June 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/helper_functions/')
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/')
ccc;

%% Set up directories and subj info

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
path_topup_fmparams = '/projectnb/somerslab/tom/projects/spacetime_network/data/fm_acqparams.txt';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.('spacetimeRuns'),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir = [projectDir, 'data/'];

experiment_name = 'rest';
fmap_prefix = '_rs';
func_data_dir = 'rest';

subjs_converted_fmaps = {'NM', 'RR', 'RT', 'NS', 'PL', 'PT', 'TP', 'UV', 'LA', 'LN', 'MK',...
    'KQ', 'SL', 'PQ'}; % these subjects had their fmaps registered to their functional data in order to convert from 2.2 - 2.3mm and will need to be handled differently

%% Loop through subjs
parfor ss = 1:length(subjCodes)

    % Set up subj varibles
    subjCode = subjCodes{ss};
    if ismember(subjCode, {'MM', 'PP', 'AH'}) % no resting state (or different sequence than normal)
        continue
    end
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    experiment_date = subjDf_cut.([experiment_name,'Date']){subjRow};
    if contains(experiment_date, '/') % different rs runs may have different dates
        experiment_date = strsplit(experiment_date, '/');
    else
        experiment_date = {experiment_date};
    end
    dirTarget = [projectDir 'data/unpacked_data_nii/' subjCode, '/' ];
    unix(['mkdir -p ' dirTarget]); % make dir if not already there

    % Check if subj with converted/registered fmaps
    if ismember(subjCode, subjs_converted_fmaps)
        reg_fmap = true;
    else
        reg_fmap = false;
    end

    % Get fieldmap run numbers
    FMruns = subjDf_cut.([experiment_name, 'FM']){subjRow};
    if contains(FMruns, '/') % different spacetime runs may have different fieldmaps
        FMruns = strsplit(FMruns, '/'); % still take all fieldmaps
    else
        FMruns = {FMruns};
    end

    for dd = 1:length(experiment_date)

        FMruns_curr = str2num(FMruns{dd});
        assert(mod(length(FMruns_curr),2)==0, ['Subj ' subjCode 'number of fieldmaps not a multiple of 2']);
        num_pairs = length(FMruns_curr)/2;
        run_inds = reshape(1:length(FMruns_curr), [2,num_pairs]);

        if reg_fmap
            num_pairs = length({dir([projectDir 'data/unpacked_data_nii/' subjCode '/rest/']).name}) - 2;
        end
        
        for pp = 1:num_pairs
            
            if reg_fmap
                fmap_runs = [pp, pp];
                fname_prefix1 = 'rs'; 
                fname_prefix2 = '_';
                fname_suffix = '.nii.gz';
            else
                fmap_runs = FMruns_curr(run_inds(:,pp));
                fname_prefix1 = 'run';
                fname_prefix2 = '_rs';
                fname_suffix = '.nii';
            end

            fmapAP_filepath = [dirTarget 'bold/sub-' subjCode '_' fname_prefix1 num2str(fmap_runs(1)) fname_prefix2 'fieldmapAP' fname_suffix];
            fmapPA_filepath = [dirTarget 'bold/sub-' subjCode '_' fname_prefix1 num2str(fmap_runs(1)) fname_prefix2 'fieldmapPA', fname_suffix];
            fmapMerged_filepath = [dirTarget 'bold/sub-' subjCode 'runs' num2str(fmap_runs(1)) num2str(fmap_runs(2)) '_rsfmapMerged.nii.gz'];
            fmaptopup_filepath = [dirTarget 'bold/sub-' subjCode 'runs' num2str(fmap_runs(1)) num2str(fmap_runs(2))  '_rsfmapTopupOut'];

            if ~isfile([fmaptopup_filepath '_fieldcoef.nii.gz'])
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
                else
                    split_outpath = split(fmapMerged_filepath,'.');
                end

                % Use topup command to convert from echo to more common fieldmaps (requires all dims of images to be even)
                unix(['topup --imain=' split_outpath{1} '.nii' ' --datain=/projectnb/somerslab/tom/projects/spacetime_network/data/fm_acqparams.txt '...
                    '--config=b02b0.cnf --out=' fmaptopup_filepath]);
            else
                disp(['Subj ' subjCode ': converted fieldmap target folder already contains this file. Remove this file if you want to re-convert. Skipping.'])
            end
        end


        %% Apply topup fieldmaps to rs functionals
        assert(mod(length(FMruns_curr),2)==0, ['Subj ' subjCode ' number of fieldmaps not a multiple of 2']);

        % Get functional scan run numbers
        func_runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
        if contains(func_runs, '/') % runs with different fieldmaps
            func_runs = replace(func_runs, '/', ','); % still take all runs
        end
        func_runs = str2num(func_runs);

        % Loop through functional runs
        for ff = 1:length(func_runs)
            if ismember(subjCode, {'AG', 'AI'})  % These subjs have rs runs 1 and 2 on one date and then 3rd rs run with on another date
                assert(num_pairs == 1, ['unexpected number of fieldmaps for subj ' subjCode ' (expected 1, got ' num2str(num_pairs)])
                topup_path = [dirTarget 'bold/sub-' subjCode 'runs' num2str(FMruns_curr(1)) num2str(FMruns_curr(2)) fmap_prefix 'fmapTopupOut'];
                if ~( (ff==3 && dd==2) || ( (ff==1 || ff==2) && dd==1 ) ) % only run this for loop if its func run 3 and date 2 or func runs 1 or 2 with date 1 bc those are the correct fieldmaps for the given func runs
                    continue
                end
            elseif reg_fmap
                assert(num_pairs == length(func_runs), ['unexpected number of fieldmaps for subj ' subjCode]);
                topup_path = [dirTarget 'bold/sub-' subjCode 'runs' num2str(ff) num2str(ff) fmap_prefix 'fmapTopupOut'];
            else
                assert(num_pairs == 1, ['unexpected number of fieldmaps for subj ' subjCode ' (expected 1, got ' num2str(num_pairs)]);
                topup_path = [dirTarget 'bold/sub-' subjCode 'runs' num2str(FMruns_curr(1)) num2str(FMruns_curr(2)) fmap_prefix 'fmapTopupOut'];
            end

            % Apply topup to funcitonal scans
            path_func = [dirTarget func_data_dir '/00' num2str(ff) '/uf.nii'];
            path_func_out = [dirTarget func_data_dir '/00' num2str(ff) '/uf_topupApplied.nii'];

            if isfile([path_func_out '.gz'])
                disp(['Subj ' subjCode ' run ' num2str(ff) ': there is already a functional file with topup applied in this subjects run directory, skipping...'])
            elseif ~isfile(path_func)
                disp(['unwarped functional (uf.nii) does not exist for subj ' subjCode ' rs run ' num2str(ff) ' ...skipping...']);
                continue
            else
                apply_topup_padding(path_func, path_topup_fmparams, topup_path, path_func_out, subjCode);
            end
        end
    end
    disp(['Finished subj ' subjCode]);
end


