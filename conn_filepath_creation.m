%%%%
% The purpose of this script is to print to console the filepaths for all
% T1s, fieldmaps, and functional runs, so that they can be copy/pasted into
% Conn Toolbox for easy intialization of data directories
%
% Created: Tom Possidente - March 2024
%%%%


addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

func_topupapplied = true; % give paths for functionals with fmaps already applied or not
resting_state = true;

if resting_state 
    data_dir = 'rest';
    experiment_name = 'rest';
else
    data_dir = 'bold';
    experiment_name = 'spacetime';
end

% Load in subject info
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;

ROI_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
struct_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/recons/';
func_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/';

% Print ROI annot paths for ROI selection
not_found = zeros(length(subjCodes),1);
for ss=1:length(subjCodes)
    
    subjROIpath = [ROI_path '/lh.' subjCodes{ss} '_all_ROIs_nomissing.annot'];
    if ~isfile(subjROIpath)
        not_found(ss) = 1;
    else
        disp(subjROIpath)
    end
end
disp(['No ROIs: ' string(subjCodes(logical(not_found)))' ])

% Change subjCodes to reflect only subjs with ROIs
if resting_state
    subjCodes = subjCodes(~logical(not_found));
    subjCodes = subjCodes(~ismember(subjCodes, {'PP','MM'})); % these two subjs have different resting state sequences
else
    subjCodes = subjCodes(~logical(not_found));
end

for ss=1:length(subjCodes)
    
    subjROIpath = [ROI_path '/' subjCodes{ss} '_all_ROIs_nomissing.surf.nii'];
    disp(subjROIpath)

end

% Print T1 paths
for ss=1:length(subjCodes)

    subjDirStruct = [struct_path '/' subjCodes{ss} '/mri/T1.nii'];
    assert(isfile(subjDirStruct), ['Subj ' subjCodes{ss} ' T1 file not found'])
    disp(subjDirStruct)

end

% Print functional paths
runs_all = {};
if func_topupapplied
    prefix = 'u';
    suffix = '_topupApplied.nii.gz';
else
    prefix = '';
    suffix = '.nii';
end
for ss=1:length(subjCodes)
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};

    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);
    runs_all{ss} = runs;

    for ii=1:length(runs)
        subjDirFunc = [func_path '/' subjCode '/' data_dir, '/00' num2str(ii) '/' prefix 'f' suffix];
        assert(isfile(subjDirFunc), ['Subj ' subjCode ' run ' num2str(runs(ii)) ' functional file not found'])
        disp(subjDirFunc)
    end

end

% Print realignment files (rp_f.txt)
for ss=1:length(subjCodes)
    subjCode = subjCodes{ss};
    for ii=1:length(runs_all{ss})
        realignment_filepath = [func_path '/' subjCode '/' data_dir '/00' num2str(ii) '/rp_f.txt'];
        assert(isfile(realignment_filepath), ['Subj ' subjCode ' run ' num2str(ii) ' realignment file not found'])
        disp(realignment_filepath)
    end
end

% Print QC files (art_regression_timeseries_auf_topupApplied.mat)
for ss=1:length(subjCodes)
    subjCode = subjCodes{ss};
    for ii=1:length(runs_all{ss})
        QC_filepath = [func_path '/' subjCode '/' data_dir '/00' num2str(ii) '/art_regression_timeseries_auf_topupApplied.mat'];
        assert(isfile(QC_filepath), ['Subj ' subjCode ' run ' num2str(ii) ' QC file not found'])
        disp(QC_filepath)
    end
end

% Print scrubbing files (art_regression_outliers_auf_topupApplied.mat)
for ss=1:length(subjCodes)
    subjCode = subjCodes{ss};
    for ii=1:length(runs_all{ss})
        scrub_filepath = [func_path '/' subjCode '/' data_dir '/00' num2str(ii) '/art_regression_outliers_auf_topupApplied.mat'];
        assert(isfile(scrub_filepath), ['Subj ' subjCode ' run ' num2str(ii) ' scrubbing file not found'])
        disp(scrub_filepath)
    end
end

% Print fieldmap paths
for ss=1:length(subjCodes)
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    experiment_date = subjDf_cut.([experiment_name,'Date']){subjRow};
    FMruns = subjDf_cut.([experiment_name, 'FM']){subjRow};

    if contains(FMruns, '/') % different spacetime runs may have different fieldmaps
        FMruns = replace(FMruns, '/', ','); % still take all fieldmaps
    end

    FMruns = str2num(FMruns);
    num_pairs = length(FMruns)/2;

    if num_pairs ~= 1 % if there is more than 1 pair of fieldmaps, things get more complicated
        if strcmp(subjCode, 'NS') % This subj has 1st run with 1st parir of FMs then next 3 runs with 2nd pair of FMs
            subjDirFM1 = [func_path '/' subjCode '/bold/sub-' subjCode 'runs' ...
                num2str(FMruns(1)) num2str(FMruns(2)) '_fmapMag.nii.gz'];
            disp(subjDirFM1)
            subjDirFM234 = [func_path '/' subjCode '/bold/sub-' subjCode 'runs' ...
                num2str(FMruns(3)) num2str(FMruns(4)) '_fmapMag.nii.gz'];
            disp(subjDirFM234)
            disp(subjDirFM234)
            disp(subjDirFM234)
        else
            error('More than 2 echo fieldmaps detected, should not happen')
        end

    else
        for ii = 1:length(runs_all{ss}) % print a fieldmap filepath for every functional run
            assert(isfile([func_path '/' subjCode '/bold/sub-' subjCode 'runs' ...
                num2str(FMruns(1)) num2str(FMruns(2)) '_fmapMag.nii.gz']), ['Subj ' subjCode ' runs ' ...
                num2str(FMruns(1)) num2str(FMruns(2)) ' Fieldmap file not found'])
            subjDirFM = [func_path '/' subjCode '/bold/sub-' subjCode 'runs' ...
                num2str(FMruns(1)) num2str(FMruns(2)) '_fmapMag.nii.gz'];
            disp(subjDirFM)
        end
    end



end
