%%% The purpose of this script is to reorganize the functional data
%%% folders into freesurfer compliant format
%%% Tom Possidente - May 2024

experiment_name = 'spacetime';


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

    %% Rename func folder
    old = [subjectsDir 'unpacked_data_nii/' num2str(subjCode) '/func'];
    new = [subjectsDir 'unpacked_data_nii/' num2str(subjCode) '/bold'];
    unix(['mv ' old ' ' new]);

    %% Rename functional folders and files
    for rr = 1:length(runs)
        func_path = [subjectsDir 'unpacked_data_nii/' num2str(subjCode) '/func/run' num2str(runs(rr)) '/'];

        % Rename run folder
        new_foldername = [subjectsDir 'unpacked_data_nii/' num2str(subjCode) '/bold/00' num2str(rr)];
        unix(['mv ' subjectsDir 'unpacked_data_nii/' num2str(subjCode) '/bold/run' num2str(runs(rr)) '/'...
              ' ' new_foldername]);
        
        % Rename functional file
        old_filename = [new_foldername '/sub-' num2str(subjCode) '_run' num2str(runs(rr)) '_spacetime.nii'];
        new_filename = [new_foldername '/f.nii'];
        unix(['mv ' old_filename ' ' new_filename]);

        % Rename event file
        old_eventfile = [old_filename(1:end-4) '_events.tsv'];
        new_eventfile = [new_foldername '/f_events.nii'];
        unix(['mv ' old_eventfile ' ' new_eventfile]);

        % Delete other files
        base = [new_foldername '/sub-' num2str(subjCode) '_run' num2str(runs(rr)) '_spacetime_'];
        unix(['rm ' base 'padded.nii']);
        unix(['rm ' base 'topupApplied.nii']);
        unix(['rm ' base 'topupApplied_events.tsv']);
        unix(['rm ' base 'topupApplied_padded.nii.gz']);
        unix(['rm ' new_foldername '/sub-' num2str(subjCode) 'run' num2str(runs(rr)) '_spacetime_topupApplied_padded.nii.gz']);

    end
end




