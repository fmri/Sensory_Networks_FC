%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create probabilistic ROIs from the
%%% passive-fixation (sensory drive) and active-passive (working memory) contrasts
%%% of subjs localizer data by finding areas of overlap
%%% Tom Possidente - August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Get Subj Codes

experiment_name = 'spacetime';

ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'RR'})); % rejected subjs
N = length(subjCodes);

%% Initialize variables
data_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/';
t_thresh = 2;
contrasts = {'ts_visual-f', 'ts_auditory-f', 'ts_tactile-f', 'ts_visual-ts_audtact', 'ts_auditory-ts_vistact', 'ts_tactile-ts_visaud'}; 
N_contrasts = length(contrasts);
hemis = {'lh', 'rh'};
N_hemis = length(hemis);
N_vertices = 163842;
prob_ROI_label = cell(N_contrasts, N_hemis);

lh_cortex_label = readtable([ROI_dir 'lh_inflated_wholecortex.label'], 'FileType','text');
rh_cortex_label = readtable([ROI_dir 'rh_inflated_wholecortex.label'], 'FileType','text');
lhrh_cortex_label = {lh_cortex_label, rh_cortex_label};

%% Loop through subjs and extract contrast tstat data

for cc = 1:N_contrasts
    contrast = contrasts{cc};
    for hh = 1:N_hemis
        all_subj_labels = [];
        for ss = 1:N
            subjCode = subjCodes{ss};

            cortex_label = lhrh_cortex_label{hh};
            hemi = hemis{hh};
            tstat_path = [data_dir subjCode '/bold/taskswitch_VAT_spacetime_contrasts_' hemi '_polyfit2hrf1/' contrast '/t.nii.gz'];
            tstat_data = MRIread(tstat_path);
            binarized_tstat_data = tstat_data.vol >= t_thresh;
            label_inds = find(binarized_tstat_data) - 1;
            label_rows = cortex_label(ismember(cortex_label.Var1, label_inds),:);
            nrows = size(all_subj_labels, 1);
            all_subj_labels(nrows+1:nrows+size(label_rows,1),:) = table2array(label_rows);
        end
        [unique_rows,~,ind] = unique(all_subj_labels,'rows');
        counts = histc(ind,unique(ind)); % count occurences of each vertex
        unique_rows(:,5) = counts;
        prob_ROI_label{cc,hh} = unique_rows;

        % Save label files with different thresholds
        for tt = 3:11
            label_fname = [ROI_dir hemi '_' contrast '_cortex_probabilistic_thresh' num2str(tt) '.label'];
            label = unique_rows(counts>=tt,:);
            label_file = fopen(label_fname,'w');
            fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label,1)) '\n']);
            writematrix(label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
            fclose(label_file);
        end

    end
end



