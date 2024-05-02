%%%%%%%%%%%%%%
% The purpose of this script is to transform the ROIs from CWB space to
% fsaverage space and split them into individual binarized, non-overlapping nii files
%%%%%%%%%%%%%%

%% Set necessary paths

vaibhav_path = '/projectnb/somerslab/vaibhavt/Projects/Trifloc/Codes';
save_out_filepath = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';

left_wbsurface = '/projectnb/somerslab/vaibhavt/HCP_Data/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.L.inflated_MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.surf.gii';
right_wbsurface = '/projectnb/somerslab/vaibhavt/HCP_Data/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedParcellation210.R.inflated_MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.surf.gii';
prob_map16 = '1wayLocalizerOnevsOthers_probabilistic_map.dscalar.nii';
nm_tstats = '1wayPilot_nm3p20_avg_tstat.dscalar.nii';
nm_wb_ROI_borders_R = 'IndividualROIs_nm3p20_rh.border';
nm_wb_ROI_borders_L = 'IndividualROIs_nm3p20_lh.border';

fsaverage_sphere164k = '/projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fsaverage_std_sphere.L.164k_fsavg_L.surf.gii';
fsaverage_sphere32k = '/projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii';


%% Disp all terminal commands necessary 
disp('load module connectomewb/1.3.2');
disp(['cd ' vaibhav_path]);

% Let's view everything - blue is visual, orange is auditory, green is tactile, red is MD
disp(['wb_view ' left_wbsurface ' ' right_wbsurface ' ' prob_map16 ' ' nm_tstats ' ' nm_wb_ROI_borders_R ...
      ' ' nm_wb_ROI_borders_L]);

% Convert border to ROI files (gi)
disp(['wb_command -border-to-rois ' left_wbsurface ' ' nm_wb_ROI_borders_L ' ' save_out_filepath 'nm_ROI_wb.L.func.gii'])
disp(['wb_command -border-to-rois ' right_wbsurface ' ' nm_wb_ROI_borders_R ' ' save_out_filepath 'nm_ROI_wb.R.func.gii'])

% Look at ROIs 
% TODO - write out all ROI label names and what modality they correspond to
% - or try to automate this process by loading in borders or wb rois file
disp(['wb_view ' left_wbsurface ' ' right_wbsurface ' ' prob_map16 ' ' nm_tstats ' ' nm_wb_ROI_borders_R ...
      ' ' nm_wb_ROI_borders_L ' ' save_out_filepath 'nm_ROI_wb.L.func.gii']);


% Convert all ROIs to fsaverage space
disp(['wb_command -metric-resample ' save_out_filepath 'nm_ROI_wb.L.func.gii ' fsaverage_sphere32k ' ' fsaverage_sphere164k ' ADAP_BARY_AREA ' save_out_filepath 'nm_ROI_fs_164.L.func.gii -area-metrics /projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /projectnb/somerslab/db/wb/standard_mesh_atlases/resample_fsaverage/fsaverage.L.midthickness_va_avg.164k_fsavg_L.shape.gii'])

% Convert gifti to nifti so freesurfer/conn can use it
disp(['mri_convert ' save_out_filepath 'nm_ROI_fs_164.L.func.gii ' save_out_filepath 'nm_ROI_fs_164.L.func.nii'])

% Now all ROIs have lost names and are not binarized

% Most likely will need to separate the nii into separate niis for each
% ROI

% Will have to delete overlap between ROIs in each subject individually

% Create table subjs by ROI (each cell has # of voxels in ROI)


