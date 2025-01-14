%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to quantify the overlap of the
%%% significant sensory drive and working memory areas for all ROIs using
%%% the localizer data from the spacetime experiment
%%%
%%% Tom Possidente - December 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load Subj Codes
reject_subjs = {'RR', 'AH', 'SL', 'AI'};
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.('spacetimeRuns'),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, reject_subjs));

%% Intialize Key Variables
data_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii_fs_localizer/';
ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
t_thresh = 2;
contrasts = {'f-vP', 'f-aP', 'f-tP' 'vA-vP', 'aA-aP', 'tA-tP'};
contrast_sets = {[1,4], [2,5], [3,6], [1,2,3], [4,5,6], [1,2,3,4,5,6]};
hemis = {'lh', 'rh'};
N = length(subjCodes);
N_contrasts = length(contrasts);
N_contrastsets = length(contrast_sets);
N_hemis = length(hemis);
N_vertices = 163842;

intersect_probs = zeros(N_vertices, N_hemis, N_contrastsets);
N_exclusives = sum(cell2mat(cellfun(@(x) length(x), contrast_sets, 'UniformOutput', false))); % each contrast in each contrast set will have an exclusive map
exclusion_probs = zeros(N_vertices, N_hemis, N_exclusives); 

lh_cortex_label = readtable([ROI_dir 'lh_inflated_wholecortex.label'], 'FileType','text');
rh_cortex_label = readtable([ROI_dir 'rh_inflated_wholecortex.label'], 'FileType','text');
lhrh_cortex_label = {lh_cortex_label, rh_cortex_label};

%% Loop through subjs
    
for ss = 1:N
    subjCode = subjCodes{ss};
    for hh = 1:N_hemis
        hemi = hemis{hh};
        contrast_data = nan(N_vertices, N_contrasts);

        % collect all contrasts tstat>thresh data
        for cc = 1:N_contrasts
            contrast = contrasts{cc};
            tstat_path = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrast '/t.nii.gz'];
            if ~isfile(tstat_path)
                disp([subjCode ' ' hemi ' ' contrast ' tstat file does not exist, skipping...' ]);
                continue
            end
            tstat_data = MRIread(tstat_path);
            if ismember(contrasts{cc}, {'f-vP', 'f-aP', 'f-tP'})
                tstat_data.vol = -tstat_data.vol; % reverse contrast for interpretability
            end
            contrast_data(:,cc) = tstat_data.vol >= t_thresh;
        end
        
        % Loop through contrast sets and get intersections/exclusions
        for ff = 1:N_contrastsets
            set = contrast_sets{ff};
            if any( all(isnan(contrast_data(:,set)),1) ) % if any of the contrasts in this set do not exist, skip this contrast set
                continue
            end
            intersect_probs(:,hh,ff) = intersect_probs(:,hh,ff) + ( sum(contrast_data(:,set),2)==length(set) );
            exclusions = find_exclusive_regions(contrast_data(:,set));
        end

        % cortex_label = lhrh_cortex_label{hh};
        % label_inds = find(binarized_tstat_data) - 1;
        % label_rows = cortex_label(ismember(cortex_label.Var1, label_inds),:);
        % nrows = size(all_subj_labels, 1);
        % all_subj_labels(nrows+1:nrows+size(label_rows,1),:) = table2array(label_rows);

    end
end



% [unique_rows,~,ind] = unique(all_subj_labels,'rows');
% counts = histc(ind,unique(ind)); % count occurences of each vertex
% unique_rows(:,5) = counts;
% prob_ROI_label{cc,hh} = unique_rows;
% 
% % Save label files with different thresholds
% for tt = 3:max(counts)-1
%     label_fname = [ROI_dir hemi '_' contrast '_cortex_probabilistic_thresh' num2str(tt) '.label'];
%     label = unique_rows(counts>=tt,:);
%     label_file = fopen(label_fname,'w');
%     fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label,1)) '\n']);
%     writematrix(label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
%     fclose(label_file);
% end










