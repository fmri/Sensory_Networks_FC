%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to plot probabilistic ROIs at various
%%% thresholds with freeview
%%% Tom Possidente - August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Set Key Variables

ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
subjCode = {'probabilistic'};

ROI = 'FO';
use_fsaverage = true;
opacity = '0.3';

ROIs = ["aINS", "preSMA", "ppreCun", "dACC", ... % multisensory
    "sPCS", "iPCS", "midIFS", "aIPS", "pIPS", "DO", "LOT", "VOT",... % visual
    "tgPCS", "cIFSG", "pAud", "CO", "FO", "cmSFG"]'; % auditory
N_ROIs = length(ROIs);
color = [[0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; ... % green
    [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; ... % blue
    [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]]; % orange
ctable = table(ROIs, color(:,1), color(:,2), color(:,3), 'VariableNames', {'ROI', 'c1', 'c2', 'c3'});

%% Loop through ROIs and get filepaths
lh_ROI_paths = {};
rh_ROI_paths = {};
ROI_files = {dir(ROI_dir).name};
ROI_files = ROI_files(contains(ROI_files, 'probabilistic') & contains(ROI_files, ROI));
for rr = 1:length(ROI_files)
    fname = ROI_files{rr};
    if contains(fname, '_lh')
        lh_ROI_paths{end+1} = [ROI_dir fname];
    elseif contains(fname, '_rh')
        rh_ROI_paths{end+1} = [ROI_dir fname];
    else
        error(['no rh or lh designation in filename: ' fname])
    end

end


%% 
contrasts = {};
ss_on = false;
freeview_screenshots(subjCode, contrasts, {lh_ROI_paths}, {rh_ROI_paths}, ctable, use_fsaverage,[],[],[],[],[],[],[],opacity,[],ss_on)




