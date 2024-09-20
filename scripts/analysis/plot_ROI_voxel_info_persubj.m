%%% The purpose of this script is to plot the voxel statistics for each ROI
%%% before and after overlap deletetion as a heatmap. Not very useful tbh
%%% Tom Possidente - May 2024

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/')
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;

%% Load in voxel info matrices
ROI_voxel_info = load([projectDir, 'data/ROI_persubj_voxel_info.mat']);

%% Plotting
close all;
figure;

plot_types = {'ROI_nvoxel_orig_tbl','ROI_nvoxel_deleted_tbl','ROI_percvoxel_del_tbl','ROI_nvoxel_final_tbl'};
n_plot_types = length(plot_types);

count = 1;
for pp = 1:n_plot_types

    lh = ROI_voxel_info.([plot_types{pp} '_lh']);
    rh = ROI_voxel_info.([plot_types{pp} '_lh']);

    subplot(n_plot_types, 2, count);
    count = count + 1;
    heatmap(lh.Properties.VariableNames(1:end-1), lh{:,end}, table2array(lh(:,1:end-1)), 'Colormap', jet);
    title([plot_types{pp}, ' rh'])

    subplot(n_plot_types, 2, count);
    count = count + 1;
    heatmap(rh.Properties.VariableNames(1:end-1), rh{:,end}, table2array(rh(:,1:end-1)), 'Colormap', jet);
    title([plot_types{pp}, ' rh'])

end


