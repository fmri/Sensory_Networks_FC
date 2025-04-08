%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create probabilistic ROIs from the
%%% labeled ROIs of 21 subjs localizer data by finding areas of overlap
%%% between subjs ROIs
%%% Tom Possidente - August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Get Subj Codes

experiment_name = 'spacetime';

ROI_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'SL', 'AH', 'RR'})); % rejected subjs
N = length(subjCodes);

%% Loop through ROIs
ROIs = {'aINS_lh', 'preSMA_lh', 'ppreCun_lh', 'dACC_lh', 'aINS_rh', 'preSMA_rh', 'ppreCun_rh', 'dACC_rh' ... % multisensory
    'sPCS_lh', 'iPCS_lh', 'midIFS_lh', 'aIPS_lh', 'pIPS_lh', 'DO_lh', 'LOT_lh', 'VOT_lh', 'sPCS_rh', 'iPCS_rh', 'midIFS_rh', 'aIPS_rh', 'pIPS_rh', 'DO_rh', 'LOT_rh', 'VOT_rh' ... % visual
    'tgPCS_lh', 'cIFSG_lh', 'pAud_lh', 'CO_lh', 'FO_lh', 'cmSFG_lh', 'tgPCS_rh', 'cIFSG_rh', 'pAud_rh', 'CO_rh', 'FO_rh', 'cmSFG_rh'}'; % auditory
N_ROIs = length(ROIs);

for rr = 1:N_ROIs
    ROI = ROIs{rr};
    ROI_allsubjs = [];
    for ss = 1:N
        subjCode = subjCodes{ss};
        label_path = [ROI_dir subjCode '_' ROI '.label'];
        if isfile(label_path)
            label_data = readtable(label_path, 'FileType', 'text');
        else
            disp(['No label found for ' subjCode ' ' ROI]);
            continue;
        end
        label_data = table2array(label_data);
        vertices = size(label_data,1);
        if vertices < 100 % if ROI is smaller than 100 vertices, do not use it in probabilistic ROI creation
            disp(['Label for ' subjCode ' ' ROI ' too small to use (<100 vertices)']);
            continue;
        end
        nrows = size(ROI_allsubjs,1);
        ROI_allsubjs(nrows+1:nrows+height(label_data), :) = label_data; % add vertices to growing array
    end
    [unique_rows,~,ind] = unique(ROI_allsubjs,'rows');
    counts = histc(ind,unique(ind)); % count occurences of each vertex
    unique_rows(:,5) = counts;
    
    % Loop through all thresholds and make ROI for each
    for tt = 1:max(counts)
        label_fname = [ROI_dir 'probabilistic_' ROIs{rr} '_thresh' num2str(tt) '.label'];
        label = unique_rows(counts>=tt,:);
        if size(label,1) > 100 % if the probabilistic ROI is smaller than 100 vertices, skip it
            label_file = fopen(label_fname,'w');
            fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label,1)) '\n']);
            writematrix(label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
            fclose(label_file);
        end
    end
end





