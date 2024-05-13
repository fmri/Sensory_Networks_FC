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
    
    func_path = [subjectsDir 'unpacked_data_nii/' num2str(subjCode) '/func/'];
    files = {dir(func_path).name};
    file_mask = contains(files,["_fieldcoef.nii.gz", "movpar.txt", "topup_log"]);
    selected_filepath = files(file_mask);

    for ii = 1:length(selected_filepath)
        new_filename = replace(selected_filepath{ii}, 'fmapMerged', 'fmapTopupOut');
        unix(['mv ' func_path selected_filepath{ii} ' ' func_path new_filename])
    end
end