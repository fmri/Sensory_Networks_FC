%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to modifiy the Conn 2nd level output file
% to replace the missing ROI data with NaNs so that TFCE can be run using
% only the data for which we have robust ROIs
%
% Tom Possidente - Jan 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Intiailize key variables
N_subj = 20;
N_ROIs = 26;
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS'};

%% Load Conn 2nd level ROI data file
conn_2nd_lvl_fpath = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_localizer_task/results/secondlevel/gPPI/AllSubjects/aA(-1).vA(1)/';
load([conn_2nd_lvl_fpath 'ROI.mat'], 'ROI');

ROI_names_Conn = ROI(1).names;
ROI_names_Conn = cellfun(@(x) x(6:end-4), ROI_names_Conn, 'UniformOutput',false);


%% Load missing ROI names
load('/projectnb/somerslab/tom/projects/spacetime_network/data/missing_ROIs.mat', 'missing_ROIs');
load('/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/replacement_ROI_list.mat', 'replacement_ROIs');
missing_ROIs = [missing_ROIs'; replacement_ROIs];

%% Convert missing ROI names into subj and ROI indices
subj_inds = [];
ROI_inds = [];
for mm = 1:length(missing_ROIs)
    missing_ROI = strsplit(missing_ROIs{mm}, '_');
    if strcmp(missing_ROI{1}, 'AI') || any(strcmp(missing_ROI{2}, {'ppreCun', 'pIPS'}))
        continue;
    end
    subj_ind_mask = ismember(subjCodes, missing_ROI{1});
    assert(sum(subj_ind_mask)==1, 'number of subjects found for this missing ROI was a value other than 1');
    subj_inds(end+1) = find(subj_ind_mask);
    if strcmp(missing_ROI(3), 'lh')
        ROI_inds(end+1) = find(ismember(ROI_names_Conn, missing_ROI(2)), 1, 'first');
    else
        ROI_inds(end+1) = find(ismember(ROI_names_Conn, missing_ROI(2)), 1, 'last');
    end
end

N_missing = length(ROI_inds);

%% Loop through ROI 2nd level results and replace missing ROIs with NaN

for ii = 1:N_ROIs
    y = ROI(ii).y;
    for mm = 1:N_missing
        if ROI_inds(mm)==ii
            y(subj_inds(mm),:,:) = NaN;
        else
            y(subj_inds(mm),ROI_inds(mm),:) = NaN;
        end
    end
    ROI(ii).y = y;
end

%% Save ROI variable back to file
save([conn_2nd_lvl_fpath 'ROI.mat'], 'ROI', '-append')   
