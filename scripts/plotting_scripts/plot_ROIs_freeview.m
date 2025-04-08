%%%%
% The purpose of this script is to plot ROIs before and after overlap
% deletion on fsaverage using freeview to see which ROIs are overlapping a
% lot
% Created: Tom Possidente - May 2024
%%%%%


addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';

data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;

fsaverage_surf_lh = '/share/pkg.8/freesurfer/7.4.1_CentOS-8/install/freesurfer/subjects/fsaverage/surf/lh.inflated';
fsaverage_surf_rh = '/share/pkg.8/freesurfer/7.4.1_CentOS-8/install/freesurfer/subjects/fsaverage/surf/rh.inflated';

%% Loop through subjs and plot freeview view of large overlaps
% Load ROI overlap info
load([data_dir 'ROI_persubj_voxel_info.mat'],'ROI_perc_voxel_overlap_lh', 'ROI_perc_voxel_overlap_rh')

N = length(subjCodes);

for ss = 1:N

    subjCode = subjCodes{ss};
    lh_ROI_overlaps = ROI_perc_voxel_overlap_lh{ss};
    if ~isempty(lh_ROI_overlaps)
        ROI_names = lh_ROI_overlaps.Properties.VariableNames;
        lh_ROI_overlaps_binarized = lh_ROI_overlaps > 0.25;
        lh_ROI_overlaps_binarized = lh_ROI_overlaps_binarized - diag(diag(table2array(lh_ROI_overlaps_binarized))); % Zero out diagonal
        [inds1, inds2] = find(table2array(lh_ROI_overlaps_binarized));
        ROIs_plot = unique({ROI_names{inds1},ROI_names{inds2}});
        % Build freeview command
        command = 'freeview -f ';
        % for rr = 1:length(ROIs_plot)
        %     command = [command fsaverage_surf_lh ':overlay=' data_dir '/ROIs/' lower(subjCode) '_ROI_fs_164_' ROIs_plot{rr} '_lh.nii:opacity=0.5 '];
        % end
        for rr = 1:1
            command = [command fsaverage_surf_lh ':overlay=' data_dir '/ROIs/' lower(subjCode) '_ROI_fs_164_FO_lh.nii:opacity=0.5 ' ...
                       fsaverage_surf_lh ':overlay=' data_dir '/ROIs/' lower(subjCode) '_ROI_fs_164_Ins_mult_lh.nii:opacity=0.5'];
        end
        disp(ROIs_plot)
        disp(command)
    end
    rh_ROI_overlaps = ROI_perc_voxel_overlap_rh{ss};
    rh_ROI_overlaps_binarized = rh_ROI_overlaps > 0.25;

end





