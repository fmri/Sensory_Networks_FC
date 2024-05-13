function apply_topup_padding(filepath_func, filepath_topup_params, filepath_topupbase, filepath_func_out, subjCode)
% APPLY_TOPUP_PADDING - this function uses the FSL function "applytopup"
% with a precheck for even image dimensions, padding the images to make the
% dims even if necessary. Then unpads after the topup is applied if
% necessary
% Tom Possidente, May 2024
%
%   Inputs:
%       filepath_func: filepath to functional nii which topup will be applied to
%       filepath_topup_params: filepath to txt file containing fieldmap acquisition parameters
%       filepath_topupbase: base of filepath to files containing topup output (these filepaths
%                           should begin with filepath_topupbase and end with "_topup.log", 
%                           "_fieldcoef.nii.gz", ".nii", and "_movpar.txt"
%       filepath_func_out: filepath to functional nii output which topup has been applied to 
%       subjCode: string to indicate subj
%   Outputs:
%       None - will create output file at filepath_func_out and potentially
%              an intermediate file in the same directory with the suffix
%              "padded". 
%


images = niftiread(filepath_func); % read in functional nii
dims = size(images);
odd_dims = mod(dims,2)==1;
padded = false;

if any(odd_dims(1:3)) % if any dims are odd
    padded = true;
    disp(['Subj ' subjCode ': One or more dims of the in the functional scan is odd, padding with 0s to make all dims even. Will remove padding after fieldmap transformation.']);
    images = padarray(images, double(odd_dims), 0, 'post'); % pad odd dim with vector of zeros to make even
    info = niftiinfo(filepath_func); % modify nii info to be consistent with new dims
    info.Description = [info.Description ' - Used Matlab to pad dimensions'];
    info.ImageSize = size(images);
    info.raw.dim(2:5) = size(images);
    filepath_func_split = split(filepath_func,'.');
    filepath_func = [filepath_func_split{1} '_padded.nii']; % create new filepath 
    niftiwrite(images, filepath_func, info); % save padded nii

    filepath_func_out_orig = filepath_func_out;
    filepath_func_out_split = split(filepath_func_out,'.'); 
    filepath_func_out = [filepath_func_out_split{1} '_padded.nii']; % create filepath with _padded suffix for applytopup to indicate it is padded
end

unix(['applytopup --imain=' filepath_func ' --datain=' filepath_topup_params ' --inindex=1 --topup=' filepath_topupbase ...
    ' --out=' filepath_func_out ' --method=jac'])

% If padded, remove padding
if padded
    images = niftiread([filepath_func_out '.gz']);
    dims = size(images);
    new_dims = dims - odd_dims; % take off 1 padded dim
    images = images(1:new_dims(1), 1:new_dims(2), 1:new_dims(3), 1:new_dims(4));
    info = niftiinfo([filepath_func_out '.gz']);
    info.Description = [info.Description ' - Used Matlab to unpad dimensions'];
    info.ImageSize = size(images);
    info.raw.dim(2:5) = size(images);
    niftiwrite(images, filepath_func_out_orig, info);
end

end

