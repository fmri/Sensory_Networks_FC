%%%%
% The purpose of this script is to apply preprocessing steps to fMRI data that has already had motion
% correction (and potentially fieldmap correction) run on it
%
% Created: Tom Possidente - May 2024
%%%%%


addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';

subjectsDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
N = length(subjCodes);
fwhm = 3; %mm

subjCodes = {'AB'}

for ss = 1:N

    %% Set up subj varibles
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);

    % Make sure subjectname file exists in subj directory
    if ~isfile([subjectsDir subjCode '/subjectname'])
        unix(['echo "' subjCode '" > ' subjectsDir subjCode '/subjectname']); % create subjectname file with subjCode inside as text (necessary for preprocc-sess)
    end

    for rr = 1:length(runs)

        % Run proprocessing without MC
        if ~isfile([subjectsDir subjCode '/bold/00' num2str(rr) '/auf_topupApplied.sm' num2str(fwhm) '.nii.gz'])
            unix(['preproc-sess -s ' subjCode ' -d ' subjectsDir ' -fsd bold' ' -per-run -fwhm ' num2str(fwhm) ' -nomc -i auf_topupApplied']);
        else
           disp(['final preproc file found for subj ' subjCode ' run ' num2str(rr) '... skipping']);
        end

    end

end






% Run brain mask creation
% if ~isfile([subjectsDir subjCode '/bold/00' num2str(rr) '/fmcpr'.nii.gz'])
%     unix(['mkbrainmask-sess -s ' subjCode ' -d ' subjectsDir]);
% else
%     disp(['brainmask file found for subj ' subjCode ' run ' num2str(rr) '... skipping']);
% end
% 
% % Run spatial smoothing
% if ~isfile([subjectsDir subjCode '/bold/00' num2str(rr) '/fmcprsm' num2str(fwhm) '.nii.gz'])
%     unix(['spatialsmooth-sess -s ' subjCode ' -d ' subjectsDir ' -i fmcpr_topupApplied -o fmcprsm' ...
%         num2str(fwhm) '_topupapplied -fwhm ' num2str(fwhm)]);
% else
%     disp(['smoothed file found for subj ' subjCode ' run ' num2str(rr) '... skipping']);
% end