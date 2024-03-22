%%%%


addpath('/projectnb/somerslab/tom/helper_functions/');
ccc;

% Load in subject info
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
%subjCodes = {'AB', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'GG', 'KW', 'LA', 'LN', 'MK', 'MM', 'NS', 'PL', 'PQ', 'PT', 'SL', 'TP', 'UV'};

struct_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/recons/';
func_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii';

for ss=1:length(subjCodes)
    
    
    subjDirStruct = [struct_path '/' subjCodes{ss} '/mri/T1.nii'];
    disp(subjDirStruct)

end

for ss=1:length(subjCodes)
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    experiment_date = subjDf_cut.([experiment_name,'Date']){subjRow};
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};

    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);
    for ii=1:length(runs)
        subjDirFunc = [func_path '/' subjCodes{ss} '/func/run' num2str(runs(ii)) '/sub-' subjCodes{ss} '_run' num2str(runs(ii)) '_spacetime.nii'];
        disp(subjDirFunc)
    end

end