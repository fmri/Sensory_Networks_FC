%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to extract the gPPI beta values from the
% CONN toolbox gPPI analysis results and compute group-level statistics on
% them
% Tom Possidente - October 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Setup analysis parameters
plot_individual_betamaps = true;
save_out = false;

task = 'auditory'; % auditory or visual
if strcmp(task, 'auditory')
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_auditory_allROIs/';
    compare_conditions = [6,9]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial
    subj_rej = {'LN', 'GG', 'TP','RT', 'LA'}; %  (subj inds: 6,8,12,16,17)
    title_str = 'aS = aT';
elseif strcmp(task, 'visual')
    ROI_dataDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_spacetime_task/results/firstlevel/gPPI_visual_allROIs/';
    compare_conditions = [2,8]; % auditory spatial=6, auditory temporal=9, 2=visual temporal, 8=visual spatial
    subj_rej = {};
    title_str = 'vT-vS';
end
subjCodes = {'MM'	'PP' 'MK' 'AB' 'AD'	'LA' 'AE' 'TP' 'NM'	'AF' 'AG' 'GG' 'UV'	'PQ' 'KQ' 'LN' 'RT'	'PT' 'PL' 'NS'};
Nsubjs = length(subjCodes);
nROIs = 13*2; % 13 ROIs per hemisphere

N = length(subjCodes);
betas = nan(nROIs,nROIs,Nsubjs,2);

%% Load gPPI betas
load([ROI_dataDir 'resultsROI_Condition00' num2str(compare_conditions(1)) '.mat'], 'Z', 'names');
betas(:,:,:,1) = Z;
cond1_names = names;
load([ROI_dataDir 'resultsROI_Condition00' num2str(compare_conditions(2)) '.mat'], 'Z', 'names');
betas(:,:,:,2) = Z;
cond2_names = names;
assert(isequal(cond1_names, cond2_names))

beta_diffs = betas(:,:,:,1) - betas(:,:,:,2);
names = cellfun(@(x) x(6:end-4), cond1_names, 'UniformOutput',false);

%% Build design matrix for LME
vbias_ROIs = {'sPCS', 'iPCS', 'midIFS'};
abias_ROIs = {'tgPCS', 'cIFSG', 'cmSFG', 'CO', 'FO'};
mult_ROIs = {'aINS', 'preSMA', 'dACC'};
hemis = repelem({'lh', 'rh'},13);

data_table = table();
for ss = 1:Nsubjs
    if ismember(subjCodes{ss}, subj_rej)
        continue;
    end
    subj = ss;
    for rr1 = 1:nROIs
        ROI1 = names{rr1};
        hemi1 = hemis{rr1};
        if ismember(ROI1, abias_ROIs)
            ROI1_type = 'abias';
        elseif ismember(ROI1, vbias_ROIs)
            ROI1_type = 'vbias';
        elseif ismember(ROI1, mult_ROIs)
            ROI1_type = 'mult';
        elseif strcmp(ROI1, 'pVis')
            ROI1_type = 'pVis';
        elseif strcmp(ROI1, 'pAud')
            ROI1_type = 'pAud';
        else
            error('ROI type unknown');
        end
        for rr2 = 1:nROIs
            if rr1 ~= rr2 % don't include same ROI to itself
                ROI2 = names{rr2};
                hemi2 = hemis{rr2};
                if ismember(ROI2, abias_ROIs)
                    ROI2_type = 'abias';
                elseif ismember(ROI2, vbias_ROIs)
                    ROI2_type = 'vbias';
                elseif ismember(ROI2, mult_ROIs)
                    ROI2_type = 'mult';
                elseif strcmp(ROI2, 'pVis')
                    ROI2_type = 'pVis';
                elseif strcmp(ROI2, 'pAud')
                    ROI2_type = 'pAud';
                else
                    error('ROI type unknown');
                end
                ROItype_order = sort({ROI1_type, ROI2_type});
                connection_type = [ROItype_order{1} '<->' ROItype_order{2}];
                beta_diff = beta_diffs(rr1,rr2,ss);
                if rr1 > rr2
                    beta_diff2 = beta_diffs(rr2,rr1,ss);
                    data_table = [data_table; {mean([beta_diff, beta_diff2]), subj, hemi1, hemi2, connection_type}];
                end
            end
        end
    end
end

data_table.Properties.VariableNames = {'beta_diff', 'subject', 'hemisphere1', 'hemisphere2', 'connection_type'};
data_table.subject = categorical(data_table.subject);
data_table.hemisphere1 = categorical(data_table.hemisphere1);
data_table.hemisphere2 = categorical(data_table.hemisphere2);
data_table.connection_type = categorical(data_table.connection_type);


%% Fit LME
tic;
lme = fitlme(data_table, ['beta_diff ~ 1 + connection_type + (1 + connection_type | subject) ' ...
    ' + (1 + connection_type | hemisphere1) + (1 + connection_type | hemisphere2)'])
toc
