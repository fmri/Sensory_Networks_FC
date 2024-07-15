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

% probabilistic map
prob_map16 = [vaibhav_path '1wayLocalizerOnevsOthers_probabilistic_map.dscalar.nii']; % probabilistic maps based on activation from 16 subjs


%% Get subj codes
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
N = length(subjCodes);

%% Set ROIs
aud_ROIs = {'pSTS', 'tgPCS', 'cIFS/G', 'FO', 'CO', 'pAud', 'cmSFG'};% ignore pSTS, and cmSFG
aud_ROIs_use = {'tgPCS', 'cIFS/G', 'FO', 'CO', 'pAud'};
tac_ROIs = {'S2', 'midIns', 'S1', 'pSTS-Tac', 'ppTac', 'PMv', 'posCingSulc'};
tac_ROIs_use = {};
vis_ROIs = {'pVis', 'SPCS', 'IPCS', 'midIFS'}; % ignore preSMA-V
vis_ROIs_use = {'pVis', 'preSMA-V', 'SPCS', 'IPCS', 'midIFS'};
mult_ROIs = {'IPCS_mult', 'Lat_par_mult', 'LatPar_mult', 'SPCS_mult', 'midIFS_mult', 'cmSFG_mult', 'Ins_mult'}; % use Ins_mult and cmSFG_mult as pure MD areas? 
mult_ROIs_use = {'cmSFG_mult', 'Ins_mult'};

all_ROIs = {'pSTS', 'tgPCS', 'cIFS/G', 'FO', 'CO', 'pAud', 'cmSFG', 'S2', 'midIns',...
    'S1', 'pSTS-Tac', 'ppTac', 'PMv', 'posCingSulc', 'pVis', 'preSMA-V', 'SPCS', 'IPCS', 'midIFS',...
    'IPCS_mult', 'Lat_par_mult', 'LatPar_mult', 'SPCS_mult', 'midIFS_mult', 'cmSFG_mult', 'Ins_mult'};
all_ROIs_use = {'tgPCS', 'cIFS/G', 'FO', 'CO', 'pAud', 'pVis', 'preSMA-V', 'SPCS', ...
                'IPCS', 'midIFS', 'cmSFG_mult', 'Ins_mult'};

% Make matrices to track number of voxels in each ROI
ROI_voxels_final_lh = NaN(N,length(all_ROIs));
ROI_voxels_orig_lh = NaN(N,length(all_ROIs));
ROI_perc_deleted_lh = NaN(N,length(all_ROIs));
ROI_voxels_final_rh = NaN(N,length(all_ROIs));
ROI_voxels_orig_rh = NaN(N,length(all_ROIs));
ROI_perc_deleted_rh = NaN(N,length(all_ROIs));

% Initialize variable to look at which ROIs overlap 
ROI_perc_voxel_overlap_lh = {};
ROI_perc_voxel_overlap_rh = {};

quality_check_plotting = false;

%% Loop through subjs
for ss = 1:length(subjCodes)

    %% Setup
    subjCode = subjCodes{ss};
    tstat_map = ['1wayPilot_' lower(subjCode) '3p20_avg_tstat.dscalar.nii']; % tstat activation map for subj
    wb_ROI_borders_R = ['IndividualROIs_' lower(subjCode) '3p20_rh.border']; % border files for ROIs for subj
    wb_ROI_borders_L = ['IndividualROIs_' lower(subjCode) '3p20_lh.border'];

    cd(vaibhav_path);

    if ~isfile(wb_ROI_borders_R)
        disp([subjCode ' ROI border file not found... skipping']);
        disp(['DONE CONVERTING ROIS FOR SUBJ: ' subjCode]);
        continue
    end

    % Let's view everything - blue is visual, orange is auditory, green is tactile, red is MD
    if quality_check_plotting
        unix(['wb_view ' left_wbsurface ' ' right_wbsurface ' ' prob_map16 ' ' tstat_map ' ' wb_ROI_borders_R ...
            ' ' wb_ROI_borders_L]);
    end

    %% Extract ROI names from cwb border file
    ROI_names_lh = extract_borderfile_ROIs(wb_ROI_borders_L);
    ROI_exists_lh = ismember(all_ROIs, ROI_names_lh);
    ROI_names_rh = extract_borderfile_ROIs(wb_ROI_borders_R);
    ROI_exists_rh = ismember(all_ROIs, ROI_names_rh);

    assert(sum(ROI_exists_lh)==length(ROI_names_lh) & sum(ROI_exists_rh)==length(ROI_names_rh), 'ROIs not in list');

    %% Convert border to ROI files (gii)
    if ~isfile([save_out_filepath lower(subjCode) '_ROI_wb.L.func.gii'])
        unix(['wb_command -border-to-rois ' left_wbsurface ' ' wb_ROI_borders_L ' ' save_out_filepath lower(subjCode) '_ROI_wb.L.func.gii']);
    end
    if ~isfile([save_out_filepath lower(subjCode) '_ROI_wb.R.func.gii'])
        unix(['wb_command -border-to-rois ' right_wbsurface ' ' wb_ROI_borders_R ' ' save_out_filepath lower(subjCode) '_ROI_wb.R.func.gii']);
    end

    % Look at ROIs
    if quality_check_plotting
        unix(['wb_view ' left_wbsurface ' ' right_wbsurface ' ' prob_map16 ' ' tstat_map ' ' wb_ROI_borders_R ...
            ' ' wb_ROI_borders_L ' ' save_out_filepath lower(subjCode) '_ROI_wb.L.func.gii']);
    end

    %% Convert all ROIs to fsaverage space
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

    %% Separate each ROI into its own nii file
    ROIname_lh_santized = cell(length(ROI_names_lh),1);
    ROIname_rh_santized = cell(length(ROI_names_rh),1);
    for ff = 1:length(ROI_names_lh)
        if contains(ROI_names_lh{ff}, '/') % replace slash with underscore so fname doesn't get messed up
            ROIname_lh_santized{ff} = strrep(ROI_names_lh{ff}, '/', '_');
        else
            ROIname_lh_santized{ff} = ROI_names_lh{ff};
        end
        if ~isfile([save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_lh_santized{ff} '_lh_binarized.nii'])

            % Take single ROI
            unix(['mri_convert --frame ' num2str(ff-1) ' ' save_out_filepath lower(subjCode) '_ROI_fs_164.L.func.nii '...
                save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_lh_santized{ff} '_lh.nii']);

            % Binarize
            unix(['mri_binarize --i ' save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_lh_santized{ff} '_lh.nii' ...
                ' --o ' save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_lh_santized{ff} '_lh_binarized.nii' ...
                ' --min 0.00001']);

            % Transform from fsaverage to subj space
            %unix(['mri_surf2surf --hemi lh --srcsubject fsaverage --srcsurfval ' ...
            %    save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_lh_santized{ff} '_lh_binarized.nii' ...
            %    ' --trgsubject ' subjCode ' --trgsurfval ' save_out_filepath lower(subjCode) ...
            %    '_ROI_fs_164_' ROIname_lh_santized{ff} '_lh_binarized_subjspace.nii']);
        end

    end

    for ff = 1:length(ROI_names_rh)
        if contains(ROI_names_rh{ff}, '/') % replace slash with underscore so fname doesn't get messed up
            ROIname_rh_santized{ff} = strrep(ROI_names_rh{ff}, '/', '_');
        else
            ROIname_rh_santized{ff} = ROI_names_rh{ff};
        end
        if ~isfile([save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_rh_santized{ff} '_rh_binarized.nii'])

            % Take single ROI
            unix(['mri_convert --frame ' num2str(ff-1) ' ' save_out_filepath lower(subjCode) '_ROI_fs_164.R.func.nii '...
                save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_rh_santized{ff} '_rh.nii']);

            % Binarize
            unix(['mri_binarize --i ' save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_rh_santized{ff} '_rh.nii' ...
                ' --o ' save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_rh_santized{ff} '_rh_binarized.nii' ...
                ' --min 0.00001']);

            % Transform from fsaverage to subj space
            % unix(['mri_surf2surf --hemi rh --srcsubject fsaverage --srcsurfval ' ...
            %     save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_rh_santized{ff} '_rh_binarized.nii' ...
            %     ' --trgsubject ' subjCode ' --trgsurfval ' save_out_filepath lower(subjCode) ...
            %     '_ROI_fs_164_' ROIname_rh_santized{ff} '_rh_binarized_subjspace.nii']);
        end
    end


    %% Delete overlap between ROIs

    % First load in data from each ROI
    fsavg_hem_voxels = 163842;
    all_ROIs_data_lh = nan(fsavg_hem_voxels, length(ROI_names_lh));
    all_ROIs_data_rh = nan(fsavg_hem_voxels, length(ROI_names_rh));
    all_ROI_info_lh = {};
    all_ROI_info_rh = {};
    for ff = 1:length(ROI_names_lh)
        ROI_read_lh = MRIread([save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_lh_santized{ff} '_lh_binarized.nii']);
        all_ROIs_data_lh(:,ff) = ROI_read_lh.vol;
        all_ROI_info_lh{ff} = ROI_read_lh;
    end
    for ff = 1:length(ROI_names_rh)
        ROI_read_rh = MRIread([save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_rh_santized{ff} '_rh_binarized.nii']);
        all_ROIs_data_rh(:,ff) = ROI_read_rh.vol;
        all_ROI_info_rh{ff} = ROI_read_rh;
    end

    % Check how much each ROI overlaps
    ROI_overlap_matrix_lh = NaN(length(ROI_names_lh), length(ROI_names_lh));
    for ii = 1:length(ROI_names_lh)
        ROI_voxels_orig_lh(ss,ismember(all_ROIs,ROI_names_lh{ii})) = sum(all_ROIs_data_lh(:,ii));
        for ff = 1:length(ROI_names_lh)
            if ismember(ROI_names_lh{ii}, all_ROIs_use) && ismember(ROI_names_lh{ff}, all_ROIs_use)
                ROI_overlap_matrix_lh(ii,ff) = sum(all_ROIs_data_lh(:,ii) .* all_ROIs_data_lh(:,ff)) / sum(all_ROIs_data_lh(:,ii));
            end
        end
    end
    ROI_perc_voxel_overlap_lh{ss} = array2table(ROI_overlap_matrix_lh,'VariableNames',ROI_names_lh);
    
    ROI_overlap_matrix_rh = NaN(length(ROI_names_rh), length(ROI_names_rh));
    for ii = 1:length(ROI_names_rh)
        ROI_voxels_orig_rh(ss,ismember(all_ROIs,ROI_names_rh{ii})) = sum(all_ROIs_data_rh(:,ii));
        for ff = 1:length(ROI_names_rh)
            if ismember(ROI_names_rh{ii}, all_ROIs_use) && ismember(ROI_names_rh{ff}, all_ROIs_use)
                ROI_overlap_matrix_rh(ii,ff) = sum(all_ROIs_data_rh(:,ii) .* all_ROIs_data_rh(:,ff)) / sum(all_ROIs_data_rh(:,ii));
            end
        end
    end
    ROI_perc_voxel_overlap_rh{ss} = array2table(ROI_overlap_matrix_rh,'VariableNames',ROI_names_rh);

    % If any overlap between ROIs in any voxel, replace with 0
    %% TODO: ROIs ARE IN DIFFERENT/WRONG ORDER FOR EACH TABLE/SUBJ %%
    ROI_use_mask_lh = ismember(ROI_names_lh, all_ROIs_use);
    ROI_use_mask_rh = ismember(ROI_names_rh, all_ROIs_use);
    for vv = 1:fsavg_hem_voxels
        if sum(all_ROIs_data_lh(vv,ROI_use_mask_lh)) > 1
            all_ROIs_data_lh(vv,ROI_use_mask_lh) = zeros(1,sum(ROI_use_mask_lh));
        end
        if sum(all_ROIs_data_rh(vv,ROI_use_mask_rh)) > 1
            all_ROIs_data_rh(vv,ROI_use_mask_rh) = zeros(1,sum(ROI_use_mask_rh));
        end
    end

    % % Track how many voxels get deleted
    % ROI_voxels_final_lh(ss,ROI_exists_lh) = sum(all_ROIs_data_lh,1);
    % ROI_perc_deleted_lh(ss,ROI_exists_lh) = (ROI_voxels_orig_lh(ss,ROI_exists_lh) - ROI_voxels_final_lh(ss,ROI_exists_lh))./ROI_voxels_orig_lh(ss,ROI_exists_lh);
    % 
    % ROI_voxels_final_rh(ss,ROI_exists_rh) = sum(all_ROIs_data_rh,1);
    % ROI_perc_deleted_rh(ss,ROI_exists_rh) = (ROI_voxels_orig_rh(ss,ROI_exists_rh) - ROI_voxels_final_rh(ss,ROI_exists_rh))./ROI_voxels_orig_rh(ss,ROI_exists_rh);

    % Replace old ROI data with new ROI data for each ROI
    for ff = 1:length(ROI_names_lh)
        all_ROI_info_lh{ff}.vol = all_ROIs_data_lh(:,ff)';
        MRIwrite(all_ROI_info_lh{ff}, [save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_lh_santized{ff} '_lh_binarized_nooverlap.nii'])
    end
    for ff = 1:length(ROI_names_rh)
        all_ROI_info_rh{ff}.vol = all_ROIs_data_rh(:,ff)';
        MRIwrite(all_ROI_info_rh{ff}, [save_out_filepath lower(subjCode) '_ROI_fs_164_' ROIname_rh_santized{ff} '_rh_binarized_nooverlap.nii'])
    end

    disp(['DONE CONVERTING ROIS FOR SUBJ: ' subjCode]);
end

% Convert ROI matrices to tables
ROI_nvoxel_orig_tbl_lh = array2table(ROI_voxels_orig_lh, 'VariableNames', all_ROIs);
ROI_nvoxel_orig_tbl_lh.subjCode = subjCodes;

% ROI_percvoxel_del_tbl_lh = array2table(ROI_perc_deleted_lh, 'VariableNames', all_ROIs);
% ROI_percvoxel_del_tbl_lh.subjCode = subjCodes;
% 
% ROI_nvoxel_final_tbl_lh = array2table(ROI_voxels_final_lh, 'VariableNames', all_ROIs);
% ROI_nvoxel_final_tbl_lh.subjCode = subjCodes;
% 
% ROI_nvoxel_deleted_tbl_lh = array2table(ROI_voxels_orig_lh.*ROI_perc_deleted_lh, 'VariableNames', all_ROIs);
% ROI_nvoxel_deleted_tbl_lh.subjCode = subjCodes;

ROI_nvoxel_orig_tbl_rh = array2table(ROI_voxels_orig_rh, 'VariableNames', all_ROIs);
ROI_nvoxel_orig_tbl_rh.subjCode = subjCodes;

% ROI_percvoxel_del_tbl_rh = array2table(ROI_perc_deleted_rh, 'VariableNames', all_ROIs);
% ROI_percvoxel_del_tbl_rh.subjCode = subjCodes;

% ROI_nvoxel_final_tbl_rh = array2table(ROI_voxels_final_rh, 'VariableNames', all_ROIs);
% ROI_nvoxel_final_tbl_rh.subjCode = subjCodes;
% 
% ROI_nvoxel_deleted_tbl_rh = array2table(ROI_voxels_orig_rh.*ROI_perc_deleted_rh, 'VariableNames', all_ROIs);
% ROI_nvoxel_deleted_tbl_rh.subjCode = subjCodes;

% save('ROI_persubj_voxel_info.mat', 'ROI_nvoxel_deleted_tbl_rh','ROI_nvoxel_final_tbl_rh', 'ROI_percvoxel_del_tbl_rh',...
%     'ROI_nvoxel_orig_tbl_rh', 'ROI_nvoxel_orig_tbl_lh', 'ROI_nvoxel_deleted_tbl_lh',...
%     'ROI_nvoxel_final_tbl_lh', 'ROI_percvoxel_del_tbl_lh','ROI_perc_voxel_overlap_lh', 'ROI_perc_voxel_overlap_rh');

save('ROI_persubj_voxel_info.mat', 'ROI_nvoxel_orig_tbl_rh', 'ROI_nvoxel_orig_tbl_lh',...
     'ROI_perc_voxel_overlap_lh', 'ROI_perc_voxel_overlap_rh');



