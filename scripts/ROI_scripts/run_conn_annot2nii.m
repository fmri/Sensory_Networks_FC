%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to convert annotation files to surface
%%% nii files that CONN can use 
%%% Tom Possidente - August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;


ROI_dir = [projectDir 'data/ROIs/'];
files = {dir(ROI_dir).name};
files_cut = files(contains(files, '_avsm_ROIs3.annot') & contains(files, 'lh.'));

for ff = 1:length(files_cut)
    
    filename = conn_annot2nii_mod([ROI_dir files_cut{ff}]);
    disp(['Finished converting annot file ' files_cut{ff}]);

end


