function command = freeview_vols(subjCode, runType, runNumber, fileName, T1, funcDir, reconDir)
%FREEVIEW - displays a freeview command line command to visualize the chosen files
%
% Inputs:
%   subjCode: char - two letter subject code that must be a directory in funcDir and reconDir
%   runType: cell of chars - indicates which functional folder to use (bold, localizer, or rest)
%   fileName: cell of chars - indicates file names within runType folder
%   T1: logical - whether to include T1 scan in freeview command (OPTIONAL) default = true
%   funcDir: char - base directory for functional files (OPTIONAL). default = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/'
%   reconDir: char - base directory for T1 file (OPTIONAL). default = '/projectnb/somerslab/tom/projects/spacetime_network/data/recons/'
%
% Tom Possidente July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up Inputs
if nargin < 5 || isempty(T1)
    T1 = true;
end

if nargin < 6 || isempty(funcDir)
    funcDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/';
end

if nargin < 7 || isempty(reconDir)
    reconDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/recons/';
end

assert(length(runType) == length(runNumber) && length(runNumber) == length(fileName), 'runType, runNumber, fileName must have the same length');
assert(all(ismember(runType, {'bold' 'localizer', 'rest'})), 'runType not recognized, must be "bold", "localizer", or "rest"');


%% Start building freeview command
command = 'freeview -v ';

if T1
    T1fpath = [reconDir subjCode '/mri/T1.nii'];
    assert(isfile(T1fpath), [T1fpath ' does not exist']);
    command = [command T1fpath ' '];
end

for ii = 1:length(runType)

    fpath = [funcDir subjCode '/' runType{ii} '/00' runNumber{ii} '/' fileName{ii}];
    assert(isfile(fpath), [fpath ' does not exist']);
    command = [command fpath ' '];

end

disp(command)



