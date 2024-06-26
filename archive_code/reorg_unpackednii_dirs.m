
%% Set up subj paths and variables

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir = [projectDir, 'data/unpacked_data_nii/'];

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
        subjdir = [subjectsDir subjCode '/bold/00' num2str(rr)];
        json_filename = ['sub-' subjCode '_run' num2str(runs(rr)) '_spacetime.json'];
        cd(subjdir)
        if length(dir(subjdir)) > 4
            unix(['find . -type f -not -name "f.nii" -not -name "f.json" -not -name "' json_filename '" -print0 | xargs -0 rm --']);
            unix('rm -dr masks');
        end
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
        subjdir = [subjectsDir subjCode '/rest/00' num2str(rr)];
        json_filename = ['sub-' subjCode '_run' num2str(runs(rr)) '_spacetime.json'];
        cd(subjdir)
        if length(dir(subjdir)) > 4
            unix(['find . -type f -not -name "f.nii" -not -name "f.json" -not -name "' json_filename '" -print0 | xargs -0 rm --']);
        end
    end
end

