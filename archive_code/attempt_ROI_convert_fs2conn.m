%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to convert fs formatted nii ROI files into
% Conn formatted nii files to be used as ROIs using conn_surf_read,
% conn_surf_write, and some reshaping
%
% Tom Possidente - June 2024
%%%%%%%%%%%%%%%%%%%%%

ccc;

%% Set up filepaths
subjCode = 'MM';
ROI = 'SPCS';
ROI_filepath_lh = ['/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/' ...
                lower(subjCode) '_ROI_fs_164_' ROI '_lh_binarized_nooverlap.nii'];
ROI_filepath_rh = ['/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/' ...
                lower(subjCode) '_ROI_fs_164_' ROI '_rh_binarized_nooverlap.nii'];
test_filepath = ['/projectnb/somerslab/tom/projects/spacetime_network/data/recons/' subjCode '/label/aparc.surf.nii'];

%% 
test_connread = conn_surf_read(test_filepath);
test_MRIread = MRIread(test_filepath);

%% Read in lh and rh ROI
conn_freesurfer_read_surf
lh_SPCS = conn_surf_read(ROI_filepath_lh);
rh_SPCS = conn_surf_read(ROI_filepath_rh);
SPCS = [lh_SPCS, rh_SPCS];
SPCS_reshape = reshape(SPCS, 83, 42, 94);
conn_surf_write('test.surf.nii', SPCS_reshape); % Now test if conn will load this as an ROI in the GUI
