%%%%
% The purpose of this script is to print to console the filepaths for all
% T1s, fieldmaps, and functional runs, so that they can be copy/pasted into
% Conn Toolbox for easy intialization of data directories
%
% Created: Tom Possidente - March 2024
%%%%


addpath('/projectnb/somerslab/tom/helper_functions/');
ccc;

% Load in subject info
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;

struct_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/recons/';
func_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii';

% Print T1 paths
for ss=1:length(subjCodes)

    subjDirStruct = [struct_path '/' subjCodes{ss} '/mri/T1.nii'];
    assert(isfile(subjDirStruct), ['Subj ' subjCodes{ss} ' T1 file not found'])
    disp(subjDirStruct)

end

% Print functional paths
runs_all = {};
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
        subjDirFunc = [func_path '/' subjCode '/func/run' num2str(runs(ii)) '/sub-' subjCode '_run' num2str(runs(ii)) '_spacetime.nii'];
        assert(isfile(subjDirFunc), ['Subj ' subjCode ' run ' num2str(runs(ii)) ' functional file not found'])
        disp(subjDirFunc)
    end

end


% Print lh.aparc.annot paths for ROI selection
for ss=1:length(subjCodes)

    subjDirStruct = [struct_path '/' subjCodes{ss} '/label/lh.aparc.annot'];
    assert(isfile(subjDirStruct), ['Subj ' subjCodes{ss} ' lh.aparc.annot file not found'])
    disp(subjDirStruct)

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
            subjDirFM1 = [func_path '/' subjCode '/func/sub-' subjCode 'runs' ...
                num2str(FMruns(1)) num2str(FMruns(2)) '_fmapMag.nii.gz'];
            disp(subjDirFM1)
            subjDirFM234 = [func_path '/' subjCode '/func/sub-' subjCode 'runs' ...
                num2str(FMruns(3)) num2str(FMruns(4)) '_fmapMag.nii.gz'];
            disp(subjDirFM234)
            disp(subjDirFM234)
            disp(subjDirFM234)
        else
            error('More than 2 echo fieldmaps detected, should not happen')
        end

    else
        for ii = 1:length(runs_all{ss}) % print a fieldmap filepath for every functional run
            assert(isfile([func_path '/' subjCode '/func/sub-' subjCode 'runs' ...
                num2str(FMruns(1)) num2str(FMruns(2)) '_fmapMag.nii.gz']), ['Subj ' subjCode ' runs ' ...
                num2str(FMruns(1)) num2str(FMruns(2)) ' Fieldmap file not found'])
            subjDirFM = [func_path '/' subjCode '/func/sub-' subjCode 'runs' ...
                num2str(FMruns(1)) num2str(FMruns(2)) '_fmapMag.nii.gz'];
            disp(subjDirFM)
        end
    end



end
