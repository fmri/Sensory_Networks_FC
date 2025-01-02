%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create ROIs constrained on the
%%% probabilistic ROIs separately for passive and wm data
%%% Tom Possidente - December 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load subjs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'RR', 'MM', 'PP'}));
N_subj = length(subjCodes);

%% Loop through ROIs and subjs to create new ROIs
ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
contrast_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii_fs_localizer/';
ROI_list = {'sPCS', 'iPCS', 'midIFS', 'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG', 'pVis', 'pAud', 'dACC', 'preSMA', 'aINS'};
passive_contrast_list = {'f-vP', 'f-vP', 'f-vP', 'f-aP', 'f-aP', 'f-aP', 'f-aP', 'f-aP', 'f-vP', 'f-aP', 'both', 'both', 'both'};
active_contrast_list = {'vA-vP', 'vA-vP', 'vA-vP', 'aA-aP', 'aA-aP', 'aA-aP', 'aA-aP', 'aA-aP', 'aA-vP', 'aA-aP', 'both', 'both', 'both'};
hemis = {'lh', 'rh'};
N_ROIs = length(ROI_list);
for rr = 1:length(ROI_list)
    for hh = 1:2

        % Load probabilistic ROI
        ROI = ROI_list{rr};
        thresh = 5;
        prob_ROI_fpath = [ROI_dir 'probabilistic_' ROI '_' hemis{hh} '_thresh' num2str(thresh) '_dilated5.label'];
        while ~isfile(prob_ROI_fpath)
            thresh = thresh - 1;
            prob_ROI_fpath = [ROI_dir 'probabilistic_' ROI '_thresh' num2str(thresh) '_dilated5.label'];
            if thresh < 3
                error(['probabilistic ROI ' ROI ' not found: ', prob_ROI_fpath])
            end
        end
        disp(['Loading ' ROI ' using probabilistic ROI: probabilistic_' ROI '_thresh' num2str(thresh) '_dilated5.label'])
        prob_ROI = readtable(prob_ROI_fpath, 'FileType', 'text');

        % Loop through subjs
        for ss = 1:N_subj

            subjCode = subjCodes{ss};
            % load passive-fixation and active-passive contrasts
            if strcmp(passive_contrast_list{rr}, 'both')
            else
                passive_fpath = [contrast_dir subjCode '/localizer/localizer_contrasts_' hemis{hh} '/' passive_contrast_list{rr} '/t.nii.gz'];
                contrast_data = MRIread(passive_fpath); % extract tstat data from contrast
                prob_ROI_tstats = contrast_data.vol(prob_ROI.Var1+1); % get tstats for vertices in prob ROI % note that inds are off by one due to MRIread so add one
                prob_ROI_tstat_tbl = table(prob_ROI.Var1, prob_ROI_tstats', 'VariableNames', {'vertex_ind','tstat'});
                new_ROI_inds = threshold_tstat_ROI(prob_ROI_tstat_tbl, threshold, min_verts);
            end
            if ~strcmp(passive_contrast_list{rr}, 'both')
            else
                active_fpath = [contrast_dir subjCode '/localizer/localizer_conotrasts_' hemis{hh} '/' active_contrast_list{rr} '/t.nii.gz'];
            end
        end
    end

end
