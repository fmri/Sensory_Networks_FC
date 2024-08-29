%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to replace missing ROIs with the most
%%% active vertices from a probabilistic ROI
%%% Tom Possidente - August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load ROI size table
data_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/';
load([data_dir 'ROI_size_table.mat'], 'ROI_size_table');
subjCodes = ROI_size_table.Row;
ROI_names = ROI_size_table.Properties.VariableNames;
ROI_names = replace(ROI_names,' ','_');
[N,N_ROIs] = size(ROI_size_table);

%% Loop through subjs and ROIs and replace missing ones based on probabilistic ROIs
ROI_dir = [data_dir 'ROIs/'];
ROI_files = {dir(ROI_dir).name};
prob_ROIs = ROI_files(contains(ROI_files,'probabilistic_'));
VA_contrast_ROIs = {'sPCS_lh', 'sPCS_rh', 'iPCS_lh', 'iPCS_rh', 'midIFS_lh', 'midIFS_rh', 'tgPCS_lh', 'tgPCS_rh', ...
    'cIFSG_lh', 'cIFSG_rh', 'CO_lh', 'CO_rh', 'FO_lh', 'FO_rh', 'cmSFG_lh', 'cmSFG_rh', 'pIPS_rh'}; % V-A contrast ROIs with more than 0 ROIs missing
AP_contrast_ROIs = {'dACC_lh', 'dACC_rh', 'ppreCun_lh', 'ppreCun_rh', 'aINS_lh'}; % A-P contrast ROIs with more than 0 ROIs missing
visual_ROIs = {'sPCS_lh', 'sPCS_rh', 'iPCS_lh', 'iPCS_rh', 'midIFS_lh', 'midIFS_rh', 'pIPS_rh'};
auditory_ROIs = {'tgPCS_lh', 'tgPCS_rh', 'cIFSG_lh', 'cIFSG_rh', 'CO_lh', 'CO_rh', 'FO_lh', 'FO_rh', 'cmSFG_lh', 'cmSFG_rh'};
mult_ROIs = AP_contrast_ROIs;

ROI_badverts = {};
ROI_num_badverts = [];
count = 0;
for ss = 1:N % loop subjs
    subjCode = subjCodes{ss};
    for rr = 1:N_ROIs % loop ROIs
        ROI_size = ROI_size_table{ss,rr};
        if isnan(ROI_size) || ROI_size<100
            count = count + 1;
            % load probabilistic ROI
            ROI = ROI_names{rr};
            thresh = 5; % N=5 threshold on probabilistic ROIs
            prob_ROI_fpath = [ROI_dir 'probabilistic_' ROI '_thresh' num2str(thresh) '_dilated5.label'];
            while ~isfile(prob_ROI_fpath)
                thresh = thresh - 1;
                prob_ROI_fpath = [ROI_dir 'probabilistic_' ROI '_thresh' num2str(thresh) '_dilated5.label'];
                if thresh < 3
                    error(['probabilistic ROI for ' subjCode ' ' ROI ' not found'])
                end
            end
            disp(['Replacing ' subjCode ' ' ROI ' using probabilistic ROI: probabilistic_' ROI '_thresh' num2str(thresh) '_dilated5.label'])
            prob_ROI = readtable(prob_ROI_fpath, 'FileType', 'text');

            if any(ismember([0,163842], prob_ROI.Var1))
                keyboard; % checking to see if label files start indices at 0
            end

            % Load in appropriate contrast
            if ismember(ROI, VA_contrast_ROIs)
                contrast = 'V-A';
            elseif ismember(ROI, AP_contrast_ROIs)
                contrast = 'A-P';
            else
                error(['Unexpected missing ROI ' ROI ' (this ROI should not be missing)']);
            end
            hemi = ROI(end-1:end);
            contrast_fpath = [data_dir 'unpacked_data_nii_fs_localizer/' subjCode '/localizer/localizer_contrasts_' hemi '/' contrast '/t.nii.gz'];
            contrast_data = MRIread(contrast_fpath); % extract tstat data from contrast

            % Extract and sort tstats to get indices for new ROI
            prob_ROI_tstats = contrast_data.vol(prob_ROI.Var1+1); % get tstats for vertices in prob ROI % note that inds are off by one due to MRIread so add one 
            prob_ROI_tstat_tbl = table(prob_ROI.Var1, prob_ROI_tstats', 'VariableNames', {'vertex_ind','tstat'});
            sorted_tstat_tbl = sortrows(prob_ROI_tstat_tbl, 'tstat', 'descend');

            if ismember(ROI,visual_ROIs) || ismember(ROI,mult_ROIs)
                top_100_tbl = sorted_tstat_tbl(1:100,:);
                if ~all(top_100_tbl.tstat>0)
                    ROI_badverts{end+1} = [subjCode '_' ROI];
                    ROI_num_badverts(end+1) = sum(top_100_tbl.tstat<0);
                end
            elseif(ismember(ROI,auditory_ROIs))
                top_100_tbl = sorted_tstat_tbl(end-99:end,:);
                if ~all(top_100_tbl.tstat<0)
                    ROI_badverts{end+1} = [subjCode '_' ROI];
                    ROI_num_badverts(end+1) = sum(top_100_tbl.tstat>0);
                end
            else
                error(['Unexpected missing ROI ' ROI ' (this ROI should not be missing)']);
            end
            
            new_ROI = table2array(prob_ROI(ismember(prob_ROI.Var1, top_100_tbl.vertex_ind),:));
            new_ROI_label_fpath = [ROI_dir subjCode '_' ROI '_replacement.label'];
            label_file = fopen(new_ROI_label_fpath,'w');
            fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(new_ROI,1)) '\n']);
            writematrix(new_ROI, new_ROI_label_fpath, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
            fclose(label_file);

        end
    end
end









