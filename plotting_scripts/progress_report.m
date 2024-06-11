%%% The purpose of this script is to plot progress reports for each subj to
%%% keep track of how far along their data processing is
%%% Tom Possidente - May 2024

%% Set up paths and variables
addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
N = length(subjCodes);

base_dir = [projectDir, 'data/'];

steps = {'copy DICOMs', 'unpack T1/task/fmaps', 'unpacked rs', 'event files', 'recon T1s', 'MC task', 'MC rs' ...
         'fieldmap topup', 'preprocessed task', 'preprocessed rs', 'ROIs made', 'Conn denoised'};

steps_per_subj = zeros(N, length(steps));

DICOMs_path = [base_dir 'copied_DICOMs/'];
dicom_dirs = {dir(DICOMs_path).name};


for ss = 1:N
    
    subjCode = subjCodes{ss};
    
    % Check for copied dicoms
    if any(contains(dicom_dirs, subjCode))
        steps_per_subj(ss,1) = 1;
    end

    % Check for unpacked t1, task functional, and fieldmaps
    hasrun1 = isfile([base_dir, 'unpacked_data_nii/' subjCode '/bold/001/f.nii']);
    dirs = {dir([base_dir, 'unpacked_data_nii/' subjCode '/anat/']).name};
    hasanat = length(dirs)==3;
    dirs = {dir([base_dir, 'unpacked_data_nii/' subjCode '/bold/']).name};
    has_fm1 = any(contains(dirs, 'fieldmapAP.nii'));
    has_fm2 = any(contains(dirs, 'fieldmapPA.nii'));
    if hasrun1 && hasanat && has_fm1 && has_fm2
        steps_per_subj(ss,2) = 1;
    end

    % Check for unpacked resting state
    has_rs1 = isfile([base_dir, 'unpacked_data_nii/' subjCode '/rest/001/f.nii']);
    if has_rs1
        steps_per_subj(ss,3) = 1;
    end

    % Check for event file
    dirs = {dir([base_dir, 'unpacked_data_nii/' subjCode '/bold/001/']).name};
    if any(contains(dirs, '_events.tsv'))
        steps_per_subj(ss,4) = 1;
    end

    % Check for T1 recon
    if isfile([base_dir, 'recons/' subjCode '/mri/T1.nii'])
        steps_per_subj(ss,5) = 1;
    end

    % Check for motion corrected task functionals
    if isfile([base_dir, 'unpacked_data_nii/' subjCode '/bold/001/fmcpr.nii.gz'])
        steps_per_subj(ss,6) = 1;
    end

    % Check for motion corrected rs functionals
    if isfile([base_dir, 'unpacked_data_nii/' subjCode '/rest/001/fmcpr.nii.gz'])
        steps_per_subj(ss,7) = 1;
    end

    % Check for topup fieldmap creation
    dirs = {dir([base_dir, 'unpacked_data_nii/' subjCode '/bold/001/']).name};
    if any(contains(dirs, 'topupApplied') & contains(dirs, '.nii'))
        steps_per_subj(ss,8) = 1;
    end

    % Check for freesurfer preprocessed task
    dirs = {dir([base_dir, 'unpacked_data_nii/' subjCode '/bold/001/']).name};
    if any(contains(dirs, 'fmcpr_topupApplied.sm'))
        steps_per_subj(ss,9) = 1;
    end

    % Check for freesurfer preprocessed rs
    dirs = {dir([base_dir, 'unpacked_data_nii/' subjCode '/rest/001/']).name};
    if any(contains(dirs, 'fmcpr_topupApplied.sm'))
        steps_per_subj(ss,10) = 1;
    end

    % Check for ROIs
    dirs = {dir([base_dir 'ROIs/']).name};
    if sum(contains(dirs, [lower(subjCode) '_ROI_fs_164_'])) > 30
        steps_per_subj(ss,11) = 1;
    end

    
end

%% Plotting progress report
close all;
figure;
n_steps = length(steps);

for ss = 1:N

    subplot(N,1,ss);
    xline(1:n_steps);
    for step = 1:n_steps
        if steps_per_subj(ss,step) == 1
            fill([step-1 step step, step-1], [0 0 1 1], 'g');
            hold on;
        end
    end
    xlim([0,n_steps]);
    ylabel(subjCodes{ss});
    xticklabels('')
    yticklabels('')
end

steps_xticklabels = [' ', steps];
xticks(0:length(steps));
xticklabels(steps_xticklabels);

