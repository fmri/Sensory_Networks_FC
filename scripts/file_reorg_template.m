experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
%subjectsDir_src = [projectDir, 'data/unpacked_data_nii/'];
subjectsDir_trg = [projectDir, 'data/unpacked_data_nii_fs_localizer/'];

%% Loop over subjs and reorg dirs
for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};

    subjDir = [subjectsDir_trg subjCode '/localizer/'];

    hem = {'lh', 'rh'};
    for hh = 1:2
        conDir = [subjDir 'localizer_contrasts_' hem{hh} '/'];
        cont_dir_contents = dir(conDir);
        subDirs = cont_dir_contents([cont_dir_contents.isdir]);
        subDirs = {subDirs(3:end).name};
        subDirs = subDirs(~ismember(subDirs, 'res'));
        for sd = 1:length(subDirs)
            unix(['cp ' conDir subDirs{sd} '/sig.nii.gz ' conDir subDirs{sd} '/' subDirs{sd} '_sig.nii.gz'])
        end
    end
    
    %subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    % runs = subjDf_cut.('x1WayLocalizerRuns'){subjRow};
    % if isempty(runs)
    %     runs = subjDf_cut.('x3WayLocalizerRuns'){subjRow};
    % end
    % if contains(runs, '/') % runs with different fieldmaps
    %     runs = replace(runs, '/', ','); % still take all runs
    % end
    % runs = str2num(runs);
    %unix(['mkdir ' subjectsDir_trg subjCode '/localizer/']);
    

    % for rr = 1:length(runs)
    %     %%%%%%%%%%%%%%%%%
    %     %unix(['mkdir ' subjectsDir_trg subjCode '/localizer/00' num2str(rr) '/']);
    %     runDir_src1 = [subjectsDir_src subjCode '/localizer/00' num2str(rr) '/f.nii'];
    %     runDir_trg1 = [subjectsDir_trg subjCode '/localizer/00' num2str(rr) '/f_events.tsv'];
    %     unix(['cp -r ' runDir_src1 ' ' runDir_trg1]);
    %     %%%%%%%%%%%%%%%
    % end
end


%% ROIs
% ROI_dir = [projectDir 'data/ROIs/'];
% files = {dir(ROI_dir).name};
% files_cut = files(contains(files, 'all_ROIs.surf.'));
% 
% for ff = 1:length(files_cut)
% 
%     disp(files_cut(ff))
%     unix(['rm ' ROI_dir files_cut{ff}]);
% 
% end

