%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to take fully preprocessed surface
% functional files from Conn format and convert them to gifti so they can
% be used in freesurfer
%
% Tom Possidente July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccc;

%%
filename = 'dsauf_topupApplied.surf.nii';
data=conn_surf_read(filename);
data=reshape(data,163842,2,[]);
save(gifti(permute(data(:,1,:),[1,3,2])),conn_prepend('lh.',filename,'.gii')); 
save(gifti(permute(data(:,2,:),[1,3,2])),conn_prepend('rh.',filename,'.gii'));
