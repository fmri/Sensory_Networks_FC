
%% Set up subj paths and variables

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir_src = [projectDir, 'data/unpacked_data_nii/'];
subjectsDir_trg = [projectDir, 'data/unpacked_data_nii_backup/'];

%% Loop over subjs and reorg dirs

for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);

    for rr = 1:length(runs)
        subjdir_src = [subjectsDir_src subjCode '/bold/00' num2str(rr) '/'];
        subjdir_trg = [subjectsDir_trg subjCode '/bold/00' num2str(rr) '/'];
        json_filename = ['sub-' subjCode '_run' num2str(runs(rr)) '_spacetime.json'];
        unix(['cp ' subjdir_src '/f.nii ' subjdir_trg 'f.nii']);
        unix(['cp ' subjdir_src json_filename ' ' subjdir_trg json_filename]);
        unix(['cp ' subjdir_src '/f.nii ' subjdir_trg 'f.json']);
        unix(['cp ' subjdir_src '/f_events.tsv' subjdir_trg 'f_events.tsv']);
    end
end



for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    runs = subjDf_cut.('restRuns'){subjRow};
    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);

    for rr = 1:length(runs)
        subjdir_src = [subjectsDir_src subjCode '/rest/00' num2str(rr) '/'];
        subjdir_trg = [subjectsDir_trg subjCode '/rest/00' num2str(rr) '/'];
        json_filename = ['sub-' subjCode '_run' num2str(runs(rr)) '_spacetime.json'];
        unix(['cp ' subjdir_src '/f.nii ' subjdir_trg 'f.nii']);
        unix(['cp ' subjdir_src json_filename ' ' subjdir_trg json_filename]);
        unix(['cp ' subjdir_src '/f.nii ' subjdir_trg 'f.json']);
    end
end

%%
% Copy over fieldmaps
for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};
    unix(['cp ' subjectsDir_src subjCode '/bold/*fieldmapAP* ' subjectsDir_trg subjCode '/bold/']);
    unix(['cp ' subjectsDir_src subjCode '/bold/*fieldmapPA* ' subjectsDir_trg subjCode '/bold/']);
    unix(['cp ' subjectsDir_src subjCode '/bold/*fmapTopupOut* ' subjectsDir_trg subjCode '/bold/']);
end




