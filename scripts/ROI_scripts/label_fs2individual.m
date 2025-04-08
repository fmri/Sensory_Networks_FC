%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to transform the freesurfer ROI labels
%%% for all subjects from fsaverage space into individual subject space to
%%% make sure they match well to the original contrasts in individual space
%%% Tom Possidente - August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Get subj IDs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
N = length(subjCodes);

%% Set up paths
base_path = '/projectnb/somerslab/tom/projects/sensory_networks_FC/';
label_path = [base_path 'data/ROIs/'];


%% Loop through subjs and labels and convert labels
label_files = {dir(label_path).name};

parfor ss = 1:N

    subjID = subjCodes{ss};
    subj_label_fname = label_files(contains(label_files, subjID)); % get all ROI labels for this subj

    for ii = 1:length(subj_label_fname)
        if contains(subj_label_fname{ii}, 'aINS') && strcmp(subjID, 'NS') && ~contains(subj_label_fname{ii}, 'NS_aINS_')
            continue
        end
        new_label_fname = [subj_label_fname{ii}(1:end-6) '_indiv.label'];
        hemisphere = subj_label_fname{ii}(end-7:end-6);
        unix(['mri_label2label --srclabel ' label_path subj_label_fname{ii} ' --srcsubject fsaverage --trglabel ' ...
            label_path new_label_fname ' --trgsubject ' subjID ' --regmethod surface --hemi ' hemisphere]) % actually transform from fsaverage to subj space
    end

    disp(['Finished converting ROI labels for ' subjID]);

end









