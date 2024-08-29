experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;


ROI_dir = [projectDir 'data/ROIs/'];
files = {dir(ROI_dir).name};
files_cut = files(contains(files, 'ROIs.annot') & contains(files, 'lh.'));

for ff = 1:length(files_cut)
    
    filename = conn_annot2nii_mod([ROI_dir files_cut{ff}]);
    disp(['Finished converting annot file ' files_cut{ff}]);

end


