%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to find out which ROIs are replacement ROIs
% and get a list of them for analysis of performance of replacement vs.
% normal ROIs
% Tom Possidente - December 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/'));
ccc;

%% Search directory for replacement ROIs
ROI_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs/';
files = {dir(ROI_dir).name};
files_cut = files(contains(files, 'replacement.label'));

replacement_ROIs = cell(length(files_cut), 1);

for ff = 1:length(files_cut)
    file = files_cut{ff};
    replacement_ROIs{ff} = file(1:end-18);
end

%% Save
save('replacement_ROI_list.mat', 'replacement_ROIs');