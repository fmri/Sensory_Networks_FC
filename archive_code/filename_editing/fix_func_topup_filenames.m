%%% Fix filenames for fmaps produced by topup. Mislabeled as
%%% '*_fmapMagPhase.nii', should actually be '*_fmapTopup.nii'
%%% Tom Possidente - May 2024

experiment_name = 'spacetime';
unpack_t1s = false;
unpack_func = false;
unpack_fieldmaps = true;
convert_fieldmaps = true; % convert fieldmaps from spin-echo to mag/phase (more commonly used)

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
dicomsBase=[projectDir 'data/copied_DICOMs/'];
path_topup_fmparams = '/projectnb/somerslab/tom/projects/spacetime_network/data/fm_acqparams.txt';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir = [projectDir, 'data/'];

%% Loop through subjs
for ss = 1:length(subjCodes)

    %% Set up subj varibles
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);


    for rr = 1:length(runs)
        func_path = [subjectsDir 'unpacked_data_nii/' num2str(subjCode) '/func/run' num2str(runs(rr)) '/'];
        files = {dir(func_path).name};
        file_mask = contains(files,'topupApplied');
        selected_filepath = files(file_mask);

        for ii = 1:length(selected_filepath)
            new_filename = replace(selected_filepath{ii}, 'run', '_run');
            unix(['mv ' func_path selected_filepath{ii} ' ' func_path new_filename])
        end
    end
end