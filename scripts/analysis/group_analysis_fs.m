%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to perform a group level activity
%%% analysis on the spacetime localizer data sensory drive and WM contrasts
%%% for visual and auditory stimulation/tasks.
%%%
%%% Tom Possidente - December 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% load subjIDs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'RR', 'MM', 'PP', 'AH'}));
N_subjs = length(subjCodes);

%% Loop through subjs and contrasts and construct mri_concat command
contrasts = {'f-vP', 'f-aP', 'vA-vP', 'aA-aP'};
N_contrasts = length(contrasts);
subj_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii_fs_localizer/';
path_suffix = '/localizer/localizer_contrasts_';
hemis = {'lh', 'rh'};

for hh = 1:length(hemis)
    hemi = hemis{hh};
    for cc = 1:N_contrasts
        cmd = 'mri_concat ';
        contrast = contrasts{cc};
        for ss = 1:N_subjs
            subjCode = subjCodes{ss};
            path = [subj_path subjCode path_suffix hemi '/' contrast '/ces.nii.gz'];
            cmd = [cmd path ' '];
        end
        cmd = [cmd '--o ' hemi '.ces.localizer_groupavg_' contrast '.nii'];
        unix(cmd);
    end
end

%% Loop through concatenated data and run glm for each contrast and hemisphere
groupavg_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii_fs_localizer/grouplevel/';

for cc = 1:N_contrasts
    contrast = contrasts{cc};
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        outdir = ['/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii_fs_localizer/grouplevel/' hemi '.ces.localizer_groupavg_' contrast '.glmres'];
        unix(['mri_glmfit --y ' groupavg_dir hemi '.ces.localizer_groupavg_' contrast '.nii --surf fsaverage ' hemi ' --osgm --glmdir ' outdir ' --cortex']);
    end
end
