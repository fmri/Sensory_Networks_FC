%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to find the distance (mm) between ROIs
% Tom Possidente - Nov 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/'));
ccc;


%% Initialize Key Variables
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
ROI_names = {'tgPCS_lh', 'FO_lh', 'CO_lh', 'cIFSG_lh', 'cmSFG_lh', 'pAud_lh', ...
    'sPCS_lh', 'iPCS_lh', 'midIFS_lh', 'pVis_lh', ...
    'sm_aINS_lh', 'sm_preSMA_lh', 'sm_dACC_lh', 'sm_sPCS_lh', 'sm_iPCS_lh', 'sm_midFSG_lh',...
    'tgPCS_rh', 'FO_rh', 'CO_rh', 'cIFSG_rh', 'cmSFG_rh', 'pAud_rh', ...
    'sPCS_rh', 'iPCS_rh', 'midIFS_rh', 'pVis_rh', ...
    'sm_aINS_rh', 'sm_preSMA_rh', 'sm_dACC_rh', 'sm_sPCS_rh', 'sm_iPCS_rh', 'sm_midFSG_rh'};
hemis = [repmat({'lh'},1,16), repmat({'rh'},1,16)];

N_subjs = length(subjCodes);
N_ROIs = length(ROI_names);
ROI_final_fs_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs_final/fsaverage/';
dists = nan(N_subjs, N_ROIs, N_ROIs);

%% Loop over subjs, ROIs, hemis and calc surface areas
for ss = 1:N_subjs
    subjCode = subjCodes{ss};
    for r1 = 1:N_ROIs
        ROI1_name = ROI_names{r1};
        ROI1_hemi = hemis{r1};
        ROI1 = readtable([ROI_final_fs_dir subjCode '/' ROI1_hemi '.' ROI1_name(1:end-3) '.label'], 'FileType','text');
        [~,medoid1] = kmedoids(ROI1{:,2:4},1);
        %centroid1 = median(ROI1{:,2:4});
        for r2 = 1:N_ROIs
            if r1 ~= r2 % do not calculate distance between self
                ROI2_name = ROI_names{r2};
                ROI2_hemi = hemis{r2};
                ROI2 = readtable([ROI_final_fs_dir subjCode '/' ROI2_hemi '.' ROI2_name(1:end-3) '.label'], 'FileType','text');
                [~,medoid2] = kmedoids(ROI2{:,2:4},1);
                %centroid2 = median(ROI2{:,2:4});
                dists(ss,r1,r2) = norm(medoid1-medoid2);
            end
        end
    end
    disp(ss)
end


save('ROI_distances.mat', 'dists', 'hemis', 'subjCodes', 'ROI_names');





