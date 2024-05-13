%%%%%%%%%%%%%%
% The purpose of this script is to transform the ROIs from CWB space to
% fsaverage space and split them into individual binarized, non-overlapping nii files
%
% this module must be loaded before opening matlab and running this script: module load connectomewb/1.3.2
%
% Created 5/24 - Tom Possidente & David Beeler
%%%%%%%%%%%%%%

%% Set necessary paths
addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');

vaibhav_path = '/projectnb/somerslab/vaibhavt/Projects/Trifloc/Codes/';
save_out_filepath = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';

% workbench 32k surfaces
left_wbsurface = '/projectnb/somerslab/vaibhavt/HCP_Data/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.inflated_MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.surf.gii';
right_wbsurface = '/projectnb/somerslab/vaibhavt/HCP_Data/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.inflated_MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.surf.gii';

% fsaverage sphere surfaces (32k and 164k)
fsaverage_sphere32k_lh = '/projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii';
fsaverage_sphere164k_lh = '/projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fsaverage_std_sphere.L.164k_fsavg_L.surf.gii';
fsaverage_shape32k_lh = '/projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii';
fsaverage_shape164k_lh = '/projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fsaverage.L.midthickness_va_avg.164k_fsavg_L.shape.gii';

fsaverage_sphere32k_rh = '/projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii';
fsaverage_sphere164k_rh = '/projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fsaverage_std_sphere.R.164k_fsavg_R.surf.gii';
fsaverage_shape32k_rh = '/projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii';
fsaverage_shape164k_rh = '/projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fsaverage.R.midthickness_va_avg.164k_fsavg_R.shape.gii';

fsaverage = '$FREESURFER_HOME/subjects/fsaverage';

% probabilistic map
prob_map16 = [vaibhav_path '1wayLocalizerOnevsOthers_probabilistic_map.dscalar.nii']; % probabilistic maps based on activation from 16 subjs


%% Get subj codes
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;

%% Set ROIs
aud_ROIs = {'pSTS', 'tgPCS', 'cIFS/G', 'FO', 'CO', 'pAud', 'cmSFG'};
tac_ROIs = {'S2', 'midIns', 'S1', 'pSTS-Tac', 'ppTac', 'PMv', 'posCingSulc'};
vis_ROIs = {'pVis', 'preSMA-V', 'SPCS', 'IPCS', 'midIFS'};
mult_ROIs = {'IPCS_mult', 'Lat_par_mult', 'LatPar_mult', 'SPCS_mult', 'midIFS_mult', 'cmSFG_mult', 'Ins_mult'};

all_ROIs = {'pSTS', 'tgPCS', 'cIFS/G', 'FO', 'CO', 'pAud', 'cmSFG', 'S2', 'midIns',...
    'S1', 'pSTS-Tac', 'ppTac', 'PMv', 'posCingSulc', 'pVis', 'preSMA-V', 'SPCS', 'IPCS', 'midIFS',...
    'IPCS_mult', 'Lat_par_mult', 'LatPar_mult', 'SPCS_mult', 'midIFS_mult', 'cmSFG_mult', 'Ins_mult'};

ROI_matrix_lh = [];
ROI_matrix_rh = [];
subjs_list = {};

quality_check_plotting = false;

%% Loop through subjs
counter = 0;
for ss = 1:length(subjCodes)

    subjCode = subjCodes{ss};
    tstat_map = ['1wayPilot_' lower(subjCode) '3p20_avg_tstat.dscalar.nii']; % tstat activation map for subj
    wb_ROI_borders_R = ['IndividualROIs_' lower(subjCode) '3p20_rh.border']; % border files for ROIs for subj
    wb_ROI_borders_L = ['IndividualROIs_' lower(subjCode) '3p20_lh.border'];

    cd(vaibhav_path);

    if ~isfile(wb_ROI_borders_R)
        disp([subjCode ' ROI border file not found... skipping'])
        continue
    end

    counter = counter + 1;

    % Let's view everything - blue is visual, orange is auditory, green is tactile, red is MD
    if quality_check_plotting
        unix(['wb_view ' left_wbsurface ' ' right_wbsurface ' ' prob_map16 ' ' tstat_map ' ' wb_ROI_borders_R ...
            ' ' wb_ROI_borders_L]);
    end

    % Extract ROI names from cwb border file
    ROI_names_lh = extract_borderfile_ROIs(wb_ROI_borders_L);
    ROI_matrix_lh(counter,:) = ismember(all_ROIs, ROI_names_lh);
    ROI_names_rh = extract_borderfile_ROIs(wb_ROI_borders_R);
    ROI_matrix_rh(counter,:) = ismember(all_ROIs, ROI_names_rh);
    subj_list{counter,:} = subjCode;

    assert(sum(ROI_matrix_lh(counter,:))==length(ROI_names_lh) & sum(ROI_matrix_rh(counter,:))==length(ROI_names_rh), 'ROIs not in list');

    % Convert border to ROI files (gi)
    if ~isfile([save_out_filepath lower(subjCode) '_ROI_wb.L.func.gii'])
        unix(['wb_command -border-to-rois ' left_wbsurface ' ' wb_ROI_borders_L ' ' save_out_filepath lower(subjCode) '_ROI_wb.L.func.gii'])
    end
    if ~isfile([save_out_filepath lower(subjCode) '_ROI_wb.R.func.gii'])
        unix(['wb_command -border-to-rois ' right_wbsurface ' ' wb_ROI_borders_R ' ' save_out_filepath lower(subjCode) '_ROI_wb.R.func.gii'])
    end

    % Look at ROIs
    if quality_check_plotting
        unix(['wb_view ' left_wbsurface ' ' right_wbsurface ' ' prob_map16 ' ' tstat_map ' ' wb_ROI_borders_R ...
            ' ' wb_ROI_borders_L ' ' save_out_filepath lower(subjCode) '_ROI_wb.L.func.gii']);
    end

    % Convert all ROIs to fsaverage space
    if ~isfile([save_out_filepath lower(subjCode) '_ROI_fs_164.L.func.nii'])
        unix(['wb_command -metric-resample ' save_out_filepath lower(subjCode) '_ROI_wb.L.func.gii ' fsaverage_sphere32k_lh...
            ' ' fsaverage_sphere164k_lh ' ADAP_BARY_AREA ' save_out_filepath lower(subjCode)...
            '_ROI_fs_164.L.func.gii -area-metrics ' fsaverage_shape32k_lh ' ' fsaverage_shape164k_lh]);
            unix(['mri_convert ' save_out_filepath lower(subjCode) '_ROI_fs_164.L.func.gii ' save_out_filepath...
                  lower(subjCode) '_ROI_fs_164.L.func.nii']);

    end

    if ~isfile([save_out_filepath lower(subjCode) '_ROI_fs_164.R.func.nii'])
        unix(['wb_command -metric-resample ' save_out_filepath lower(subjCode) '_ROI_wb.R.func.gii ' fsaverage_sphere32k_rh...
            ' ' fsaverage_sphere164k_rh ' ADAP_BARY_AREA ' save_out_filepath lower(subjCode)...
            '_ROI_fs_164.R.func.gii -area-metrics ' fsaverage_shape32k_rh ' ' fsaverage_shape164k_rh]);
            unix(['mri_convert ' save_out_filepath lower(subjCode) '_ROI_fs_164.R.func.gii ' save_out_filepath ...
                  lower(subjCode) '_ROI_fs_164.R.func.nii']);
    end

    % Separate each ROI into its own nii file
    for ff = 1:length(ROI_names_lh)
        unix(['mri_convert --frame ' num2str(ff) ' ' save_out_filepath lower(subjCode) '_ROI_fs_164.L.func.nii '...
            save_out_filepath lower(subjCode) '_ROI_fs_164_' ROI_names_lh{ff} '_lh.nii']);

        % Transform from fsaverage to subj space
        unix(['mri_surf2surf --hemi lh --srcsubject ' fsaverage , ' --srcsurfval ' ...
              save_out_filepath lower(subjCode) '_ROI_fs_164_' ROI_names_lh{ff} '_lh.nii' ...
              ' --trgsubject ' subjCode ' --trgsurfval ' ...
              save_out_filepath lower(subjCode) '_ROI_fs_164_' ROI_names_lh{ff} '_lh_subjspace.nii'])
    end

    for ff = 1:length(ROI_names_rh)
        unix(['mri_convert --frame ' num2str(ff) ' ' save_out_filepath lower(subjCode) '_ROI_fs_164.R.func.nii '...
            save_out_filepath lower(subjCode) '_ROI_fs_164_' ROI_names_rh{ff} '_rh.nii']);
    end

    % % Now all ROIs have lost names and are not binarized
    %
    % % Most likely will need to separate the nii into separate niis for each
    % % ROI
    %
    % % Will have to delete overlap between ROIs in each subject individually
    %
    % % add number of vertices to ROI_table
end

% Convert ROI matrices to tables
ROI_table_lh = array2table(ROI_matrix_lh, 'VariableNames', all_ROIs);
ROI_table_lh.subjCode = subj_list;
ROI_table_rh =  array2table(ROI_matrix_rh, 'VariableNames', all_ROIs);
ROI_table_rh.subjCode = subj_list;
