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
subj_thresh = 5; % threshold for number of subjs to make probabilistic ROI
contrasts = {'f-vP', 'f-aP', 'f-tP' 'vA-vP', 'aA-aP', 'tA-tP'};
contrast_sets = {[1,4], [2,5], [3,6], [1,2,3], [4,5,6], [1,2,3,4,5,6]};

% Create label filenames for each contrast set intersection and exclusive
count = 1;
for cc = 1:length(contrast_sets)
    name_suffix = cell2mat(cellfun(@(x) [x '+'], contrasts(contrast_sets{cc}), 'UniformOutput', false));
    intersect_names{cc} = ['intersect_' name_suffix(1:end-1)];
    for ee = 1:length(contrast_sets{cc})
        contrast_set_curr = contrast_sets{cc};
        all_except_current = contrasts(contrast_set_curr(contrast_set_curr ~= contrast_set_curr(ee)));
        name_suffix = cell2mat(cellfun(@(x) ['_exclude_' x], all_except_current, 'UniformOutput', false))
        exclusive_names{count} = [contrasts{contrast_set_curr(ee)} name_suffix];
        count = count + 1;
    end
end

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
        count = 1;
        for ff = 1:N_contrastsets
            set = contrast_sets{ff};
            if any( all(isnan(contrast_data(:,set)),1) ) % if any of the contrasts in this set do not exist, skip this contrast set
                count = count + length(set);
                continue
            end
            intersect_probs(:,hh,ff) = intersect_probs(:,hh,ff) + ( sum(contrast_data(:,set),2)==length(set) );
            set_inds = count:count+length(set)-1;
            exclusion_probs(:,hh,set_inds) = squeeze(exclusion_probs(:,hh,set_inds)) + find_exclusive_regions(contrast_data(:,set));
            count = count + length(set);
        end 



    end
end

% Loop through sets and make probabilistic ROIs
for ff = 1:N_contrastsets
    set = contrast_sets{ff};
    set_inds = count:count+length(set)-1;
    for hh = 1:N_hemis

        % Make intersection label files
        intersect_prob = intersect_probs(:,hh,ff);
        cortex_label = lhrh_cortex_label{hh};
        label_inds = find(intersect_prob >= subj_thresh) - 1; % subtract 1 to make inds line up with label inds
        label = cortex_label(ismember(cortex_label.Var1, label_inds),:);
        label{:,5} = intersect_prob(intersect_prob >= subj_thresh);
        
        label_fname = [ROI_dir hemis{hh} '_' intersect_names{ff} '_probabilistic_thresh' num2str(subj_thresh) '.label'];
        label_file = fopen(label_fname,'w');
        fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label,1)) '\n']);
        writematrix(table2array(label), label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
        fclose(label_file);

        % Make exclusion label files
        for ee = 1:length(exclusive_names)
            exclusion_prob = exclusion_probs(:,hh,ee);
            label_inds = find(exclusion_prob >= subj_thresh) - 1; % subtract 1 to make inds line up with label inds
            label = cortex_label(ismember(cortex_label.Var1, label_inds),:);
            label{:,5} = exclusion_prob(exclusion_prob >= subj_thresh);
            
            label_fname = [ROI_dir hemis{hh} '_' exclusive_names{ee} '_probabilistic_thresh' num2str(subj_thresh) '.label'];
            label_file = fopen(label_fname,'w');
            fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label,1)) '\n']);
            writematrix(table2array(label), label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
            fclose(label_file);
        end

    end
end










