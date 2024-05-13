conn_filepath = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/archive_connfiles/';

load([conn_filepath '/conn_spacetime_noFmaps/results/firstlevel/SBC_01/resultsROI_Subject001_Condition001.mat']);
noFmaps = Z;

load([conn_filepath '/conn_spacetime_topupfmaps/results/firstlevel/SBC_01/resultsROI_Subject001_Condition001.mat']);
topupfmaps_connMaps = Z;

load([conn_filepath '/spacetime_conn_subj1_fieldmaps_allcondDenoise/results/firstlevel/SBC_01/resultsROI_Subject001_Condition001.mat']);
fmaps_allcondDenoise = Z;

names_short = {'inf parietal', 'lat occipital', 'mid temporal', 'pericalcarine', 'sup parietal'};

figure;
heatmap(names_short, names_short, noFmaps, 'ColorLimits', [-1,1], 'Colormap', jet);
title('MM Cond1 Conn Matrix - no fmaps')

figure;
heatmap(names_short, names_short, topupfmaps_connMaps, 'ColorLimits', [-1,1], 'Colormap', jet);
title('MM Cond1 Conn Matrix - fmaps using topup')

figure;
heatmap(names_short, names_short, fmaps_allcondDenoise, 'ColorLimits', [-1,1], 'Colormap', jet);
title('MM Cond1 Conn Matrix - fmaps using conn (denoising using all conditions')