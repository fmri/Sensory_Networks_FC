%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to copy fieldmaps that have been registered
% from 2.2 to 2.3mm isotropic by Vaibhav into my project's data structure
% Tom Possidente - June 26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccc;

%% Setup
subj_list = {'NM', 'RR', 'RT', 'NS', 'PL', 'PT', 'TP', 'UV', 'LA', 'LN', 'MK',...
    'KQ', 'SL', 'PQ'}; % these are the subjects that should have converted fieldmaps

base_dir_trg = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/';
base_dir_src = '/projectnb/somerslab/hcp_pipeline_subjects/';

for ss = 1:length(subj_list)

    subjID = subj_list{ss};
    num_rsruns = length({dir([base_dir_trg subjID '/rest/']).name}) - 2;
    for rr = 1:num_rsruns
        src_AP = [base_dir_src lower(subjID) '3p20/unprocessed/3T/rsfMRI' num2str(rr) '/SpinEchoFieldMap_AP.nii.gz'];
        src_PA = [base_dir_src lower(subjID) '3p20/unprocessed/3T/rsfMRI' num2str(rr) '/SpinEchoFieldMap_PA.nii.gz'];
        trg_AP = [base_dir_trg subjID '/bold/sub-' subjID '_rs' num2str(rr) '_fieldmapAP.nii.gz'];
        trg_PA = [base_dir_trg subjID '/bold/sub-' subjID '_rs' num2str(rr) '_fieldmapPA.nii.gz'];

        unix(['cp ' src_AP ' ' trg_AP])
        unix(['cp ' src_PA ' ' trg_PA])
    end
end