% Copy stimulus timing files to subject's directory

addpath('/projectnb/somerslab/tom/helper_functions/')
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';
unpack_t1s = true;

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
dicomsBase=[projectDir 'data/copied_DICOMs/'];

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir = [projectDir 'data/'];
stimtimeDir = [subjectsDir 'behavioral/stim_timing/'];

%%

for ss = 1:length(subjCodes)

    %% Set up subj varibles
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));

    % Get run numbers of func data
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);

    for rr = 1:length(runs)
        sourceFile = [stimtimeDir subjCode '_run' num2str(rr) '.txt'];
        destFile = [subjectsDir 'unpacked_data_nii/' subjCode '/func/run' num2str(runs(rr)) '/sub-' subjCode '_run' num2str(runs(rr)) '_spacetime_events.tsv'];
        %unix(['rm ' destFile])
        unix(['cp ' sourceFile ' ' destFile])
    end


end








