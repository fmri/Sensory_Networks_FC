function command = freeview_contrasts(subjCode, hemisphere, sig_minmax, contrast_dir, annot_fpath, fsaverage_fpath)
%freeview_contrasts - gives the freeview command to view all surface
%                     activation contrasts for a given subject
%
% Inputs:
%   subjCode: char - 2 letter subject code
%   hemisphere: char - 'lh' or 'rh' to indicate hemisphere to visualize
%   sig_minmax: numeric - minimum and maximum colorscale values for significance to visualize
%               (number indicates the number of zeros before the pvalue). Default=[1 5]
%   contrast_dir: char - directory in which to find the contrast data. Default = ['/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii_fs_localizer/' subjCode '/localizer/localizer_contrasts_' hemisphere '/']
%   annot_fpath: char - directory in which to find the annotation file for ROIs from Vaibhav. Default = ['/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/' hemisphere '.' subjCode '_all_ROIs.annot']
%   fsaverage_fpath: char - directory in which to find fsaverage inflated surface. Default = ['/projectnb/somerslab/tom/projects/spacetime_network/data/recons/fsaverage/surf/' hemisphere '.inflated']
%
% Outputs:
%   command: char - freeview command to view contrasts
%
% Tom Possidente July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/';

%% Set up Inputs
assert(nargin>=2, 'requires at least 2 input arguments (subjCode, hemisphere)');
assert(ismember(hemisphere, {'lh','rh'}), 'hemisphere argument must be "lh" or "rh"');

if nargin < 3 || isempty(sig_minmax)
    sig_minmax = [1.3 5];
end

if nargin < 4 || isempty(contrast_dir)
    %     if ismember(subjCode, {'PP'})
    %         contrast_dir = [dataDir 'unpacked_data_nii_fs_localizer/' subjCode '/localizer/localizer_contrasts_' hemisphere '_nofix/'];
    %     elseif ismember(subjCode, {'RR', 'MM'})
    %         contrast_dir = [dataDir 'unpacked_data_nii_fs_localizer/' subjCode '/localizer/localizer_contrasts_' hemisphere '_3way/'];
    %     else
    %         contrast_dir = [dataDir 'unpacked_data_nii_fs_localizer/' subjCode '/localizer/localizer_contrasts_' hemisphere '/'];
    %     end
    contrast_dir = [dataDir 'unpacked_data_nii_fs_localizer/' subjCode '/localizer/localizer_contrasts_' hemisphere '/'];

end

if nargin < 5 || isempty(annot_fpath)
    annot_fpath = [dataDir 'ROIs/' hemisphere '.' subjCode '_40_ROIs_nomissing.annot'];
end

if nargin < 6 || isempty(fsaverage_fpath)
    fsaverage_fpath = [dataDir 'recons/fsaverage/surf/' hemisphere '.inflated'];
end

%% Extract contrast surface files
cont_dir_contents = dir(contrast_dir);
subDirs = cont_dir_contents([cont_dir_contents.isdir]);
subDirs = {subDirs(3:end).name};
subDirs = subDirs(~ismember(subDirs, 'res'));

%% Start building freeview command
command = ['freeview -f ' fsaverage_fpath];

for ii = 1:length(subDirs)

    fpath = [contrast_dir subDirs{ii} '/' subDirs{ii} '_sig.nii.gz'];
    assert(isfile(fpath), [fpath ' does not exist']);
    command = [command ':overlay=' fpath ':overlay_threshold=' ...
               num2str(sig_minmax(1)) ',' num2str(sig_minmax(2))];

end

if isfile(annot_fpath)
    command = [command ':annot=' annot_fpath];
end

disp(command)



