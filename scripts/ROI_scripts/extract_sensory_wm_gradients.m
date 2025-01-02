%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to extract sensory drive and WM activity
%%% values from probabilistic ROIs for each subj to see if there is a
%%% posterior to anterior gradient.
%%% Tom Possidente - December 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load subjs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);

%% Load probabilistic ROIs
ROI_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
[~, ROI_labels{1}, ROI_info{1}] = read_annotation('/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/lh.probabilistic_ROIs.annot');
[~, ROI_labels{2}, ROI_info{2}] = read_annotation('/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/rh.probabilistic_ROIs.annot'); 
cortex_RAS{1} = readtable([ROI_path 'lh_inflated_wholecortex.label'], 'FileType', 'text');
cortex_RAS{2} = readtable([ROI_path 'rh_inflated_wholecortex.label'], 'FileType', 'text');

%% Loop through ROIs and subjs to extract data
localizer = true;
file_name = 't.nii.gz'; 
use_replacement_ROIs = false;
if localizer
    fsd = 'localizer';
    contrast_keyname = {'localizer_contrasts', ''};
    contrast_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii_fs_localizer/';
    passive_contrast_list = {'f-vP', 'f-vP', 'f-vP', 'f-aP', 'f-aP', 'f-aP', 'f-aP', 'f-aP', 'f-vP', 'f-aP', 'both', 'both', 'both'};
    passive_contrast_list_plot = {'vP-f', 'vP-f', 'vP-f', 'aP-f', 'aP-f', 'aP-f', 'aP-f', 'aP-f', 'vP-f', 'aP-f', 'both', 'both', 'both'};
    active_contrast_list = {'vA-vP', 'vA-vP', 'vA-vP', 'aA-aP', 'aA-aP', 'aA-aP', 'aA-aP', 'vA-vP', 'aA-aP', 'aA-aP', 'both', 'both', 'both'};
    subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'RR', 'MM', 'PP', 'AH'}));
    N_subj = length(subjCodes);
else
    fsd = 'bold';
    contrast_keyname = {'spacetime_contrasts', '_newcondfiles'};
    passive_contrast_list = {'pV-f', 'pV-f', 'pV-f', 'pA-f', 'pA-f', 'pA-f', 'pA-f', 'pA-f', 'pV-f', 'pA-f', 'both', 'both', 'both'};
    passive_contrast_list_plot = {'pV-f', 'pV-f', 'pV-f', 'pA-f', 'pA-f', 'pA-f', 'pA-f', 'pA-f', 'pV-f', 'pA-f', 'both', 'both', 'both'};
    active_contrast_list = {'sVtV-pV', 'sVtV-pV', 'sVtV-pV', 'sAtA-pA', 'sAtA-pA', 'sAtA-pA', 'sAtA-pA', 'sAtA-pA', 'sVtV-pV', 'sAtA-pA', 'both', 'both', 'both'};
    contrast_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/';
    subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'RR', 'AH', 'SL', 'AI'}));
    N_subj = length(subjCodes);
end
load('/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/replacement_ROI_list.mat', 'replacement_ROIs');
ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
ROI_list = {'sPCS', 'iPCS', 'midIFS', 'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG', 'pVis', 'pAud', 'dACC', 'preSMA', 'aINS'};

hemis = {'lh', 'rh'};
N_ROIs = length(ROI_list);
RAS_coords = cell(N_ROIs-3, 2);
tstats_act = cell(N_ROIs-3, 2);
tstats_pass = cell(N_ROIs-3, 2);
ROI_groupavg_tstats = cell(N_ROIs-3, 2);
RAS_coords_MD = cell(3, 2);
tstats_act_MD = cell(3, 2, 2);
tstats_pass_MD = cell(3, 2, 2);
ROI_groupavg_tstats_MD = cell(3, 2, 2);

for rr = 1:length(ROI_list)
    for hh = 1:2

        % extract probabilistic ROI
        ROI = ROI_list{rr};
        ROI_label_num = ROI_info{hh}.table(ismember(ROI_info{hh}.struct_names, ROI),5);
        vertex_inds = find(ismember(ROI_labels{hh},ROI_label_num))-1; % off by one due to vertex inds starting at 0 but indexing starting at 1
        prob_ROI = cortex_RAS{hh}(ismember(cortex_RAS{hh}{:,1},vertex_inds), :);

        if rr<=10
            RAS_coords{rr,hh} = prob_ROI;
        else
            RAS_coords_MD{rr-10,hh} = prob_ROI;
        end

        % Loop through subjs
        for ss = 1:N_subj

            subjCode = subjCodes{ss};

            % Skip if you don't want replacement ROIs
            if ~use_replacement_ROIs && ismember([subjCode '_' ROI '_' hemis{hh}], replacement_ROIs)
                continue;
            end

            % load passive-fixation and active-passive contrasts
            if strcmp(passive_contrast_list{rr}, 'both')
                assert(strcmp(active_contrast_list{rr},'both'), 'passive_contrast_list indicates MD ROI by using "both", but active_contrast_list does not');

                % Get both visual and auditory passive condition tstats
                passive_fpath = [contrast_dir subjCode '/' fsd '/' contrast_keyname{1} '_' hemis{hh} contrast_keyname{2} ...
                    '/' passive_contrast_list{1} '/' file_name];
                contrast_data = MRIread(passive_fpath); % extract tstat data from contrast
                prob_ROI_tstats = contrast_data.vol(prob_ROI.Var1+1); % get tstats for vertices in prob ROI % note that inds are off by one due to MRIread so add one
                prob_ROI_pass = table(prob_ROI.Var1, prob_ROI.Var2, prob_ROI.Var3, prob_ROI.Var4, -prob_ROI_tstats', 'VariableNames', {'vertex_ind','R','A','S','tstat'});
                tstats_pass_MD{rr-10,hh,1} = [tstats_pass_MD{rr-10,hh,1}, prob_ROI_pass.tstat];

                passive_fpath = [contrast_dir subjCode '/' fsd '/' contrast_keyname{1} '_' hemis{hh} contrast_keyname{2} ...
                    '/' passive_contrast_list{4} '/' file_name];
                contrast_data = MRIread(passive_fpath); % extract tstat data from contrast
                prob_ROI_tstats = contrast_data.vol(prob_ROI.Var1+1); % get tstats for vertices in prob ROI % note that inds are off by one due to MRIread so add one
                prob_ROI_pass = table(prob_ROI.Var1, prob_ROI.Var2, prob_ROI.Var3, prob_ROI.Var4, -prob_ROI_tstats', 'VariableNames', {'vertex_ind','R','A','S','tstat'});
                tstats_pass_MD{rr-10,hh,2} = [tstats_pass_MD{rr-10,hh,2}, prob_ROI_pass.tstat];

                % Get both visual and auditory active condition tstats
                active_fpath = [contrast_dir subjCode '/' fsd '/' contrast_keyname{1} '_' hemis{hh} contrast_keyname{2} ...
                    '/' active_contrast_list{1} '/' file_name];
                contrast_data = MRIread(active_fpath); % extract tstat data from contrast
                prob_ROI_tstats = contrast_data.vol(prob_ROI.Var1+1); % get tstats for vertices in prob ROI % note that inds are off by one due to MRIread so add one
                prob_ROI_act = table(prob_ROI.Var1, prob_ROI.Var2, prob_ROI.Var3, prob_ROI.Var4, prob_ROI_tstats', 'VariableNames', {'vertex_ind','R','A','S','tstat'});
                tstats_act_MD{rr-10,hh,1} = [tstats_act_MD{rr-10,hh,1}, prob_ROI_act.tstat];

                active_fpath = [contrast_dir subjCode '/' fsd '/' contrast_keyname{1} '_' hemis{hh} contrast_keyname{2} ...
                    '/' active_contrast_list{4} '/' file_name];
                contrast_data = MRIread(active_fpath); % extract tstat data from contrast
                prob_ROI_tstats = contrast_data.vol(prob_ROI.Var1+1); % get tstats for vertices in prob ROI % note that inds are off by one due to MRIread so add one
                prob_ROI_act = table(prob_ROI.Var1, prob_ROI.Var2, prob_ROI.Var3, prob_ROI.Var4, prob_ROI_tstats', 'VariableNames', {'vertex_ind','R','A','S','tstat'});
                tstats_act_MD{rr-10,hh,2} = [tstats_act_MD{rr-10,hh,2}, prob_ROI_act.tstat];

            else
                passive_fpath = [contrast_dir subjCode '/' fsd '/' contrast_keyname{1} '_' hemis{hh} contrast_keyname{2} ...
                    '/' passive_contrast_list{rr} '/' file_name];
                contrast_data = MRIread(passive_fpath); % extract tstat data from contrast
                prob_ROI_tstats = contrast_data.vol(prob_ROI.Var1+1); % get tstats for vertices in prob ROI % note that inds are off by one due to MRIread so add one
                prob_ROI_pass = table(prob_ROI.Var1, prob_ROI.Var2, prob_ROI.Var3, prob_ROI.Var4, -prob_ROI_tstats', 'VariableNames', {'vertex_ind','R','A','S','tstat'});
                tstats_pass{rr,hh} = [tstats_pass{rr,hh}, prob_ROI_pass.tstat];

                active_fpath = [contrast_dir subjCode '/' fsd '/' contrast_keyname{1} '_' hemis{hh} contrast_keyname{2} ...
                    '/' active_contrast_list{rr} '/' file_name];
                contrast_data = MRIread(active_fpath); % extract tstat data from contrast
                prob_ROI_tstats = contrast_data.vol(prob_ROI.Var1+1); % get tstats for vertices in prob ROI % note that inds are off by one due to MRIread so add one
                prob_ROI_act = table(prob_ROI.Var1, prob_ROI.Var2, prob_ROI.Var3, prob_ROI.Var4, prob_ROI_tstats', 'VariableNames', {'vertex_ind','R','A','S','tstat'});
                tstats_act{rr,hh} = [tstats_act{rr,hh}, prob_ROI_act.tstat];
            end
        end
        if rr<=10
            ROI_groupavg_tstats{rr,hh} = table(RAS_coords{rr,hh}.Var1, RAS_coords{rr,hh}.Var2, RAS_coords{rr,hh}.Var3, RAS_coords{rr,hh}.Var4, ...
                mean(tstats_pass{rr,hh},2), mean(tstats_act{rr,hh},2), 'VariableNames', {'index', 'R', 'A', 'S', 'tstat_pass', 'tstat_act'});
        else
            ROI_groupavg_tstats_MD{rr-10,hh,1} = table(RAS_coords_MD{rr-10,hh}.Var1, RAS_coords_MD{rr-10,hh}.Var2, RAS_coords_MD{rr-10,hh}.Var3, RAS_coords_MD{rr-10,hh}.Var4, ...
                mean(tstats_pass_MD{rr-10,hh,1},2), mean(tstats_act_MD{rr-10,hh,1},2), 'VariableNames', {'index', 'R', 'A', 'S', 'tstat_pass', 'tstat_act'});
            ROI_groupavg_tstats_MD{rr-10,hh,2} = table(RAS_coords_MD{rr-10,hh}.Var1, RAS_coords_MD{rr-10,hh}.Var2, RAS_coords_MD{rr-10,hh}.Var3, RAS_coords_MD{rr-10,hh}.Var4, ...
                mean(tstats_pass_MD{rr-10,hh,2},2), mean(tstats_act_MD{rr-10,hh,2},2), 'VariableNames', {'index', 'R', 'A', 'S', 'tstat_pass', 'tstat_act'});
        end
    end
end

%% Save out results
save('probROI_sensory_WM_tstats_localizer.mat', 'ROI_groupavg_tstats_MD', 'ROI_groupavg_tstats', 'RAS_coords_MD', 'RAS_coords', 'tstats_act', 'tstats_pass', ...
                                      'tstats_pass_MD', 'tstats_act_MD', 'ROI_list', 'subjCodes', 'hemis', 'active_contrast_list', 'passive_contrast_list')

%% Plot ROI Gradients
for rr = 1:N_ROIs
    for hh = 1:length(hemis)
        if rr<=10
            ROI_table = ROI_groupavg_tstats{rr,hh};
            min_max = [min([ROI_table.tstat_pass; ROI_table.tstat_act]), max([ROI_table.tstat_pass; ROI_table.tstat_act])];

            figure; scatter3(ROI_table,'R', 'A', 'S', 'filled', 'ColorVariable', 'tstat_pass'); colorbar;
            clim([min_max(1), min_max(2)]);
            title([ROI_list{rr} ' ' hemis{hh} ' tstat map: ' passive_contrast_list_plot{rr}]);

            figure; scatter3(ROI_table, 'R', 'A', 'S', 'filled', 'ColorVariable', 'tstat_act'); colorbar;
            clim([min_max(1), min_max(2)]);
            title([ROI_list{rr} ' ' hemis{hh} ' tstat map: ' active_contrast_list{rr}]);

            ROI_table.act_pass_tstat_diff = ROI_table.tstat_act - ROI_table.tstat_pass;
            figure; scatter3(ROI_table, 'R', 'A', 'S', 'filled', 'ColorVariable', 'act_pass_tstat_diff'); colorbar;
            title([ROI_list{rr} ' ' hemis{hh} ' tstat map: (' active_contrast_list{rr} ') - (' passive_contrast_list_plot{rr} ')']);

            keyboard;
        else
            ROI_table = ROI_groupavg_tstats_MD{rr-10,hh,1};
            min_max = [min([ROI_table.tstat_pass; ROI_table.tstat_act]), max([ROI_table.tstat_pass; ROI_table.tstat_act])];

            figure; 
            scatter3(ROI_table,'R', 'A', 'S', 'filled', 'ColorVariable', 'tstat_pass'); colorbar;
            clim([min_max(1), min_max(2)]);
            title([ROI_list{rr} ' ' hemis{hh} ' tstat map: ' passive_contrast_list_plot{1}]);

            figure; scatter3(ROI_table, 'R', 'A', 'S', 'filled', 'ColorVariable', 'tstat_act'); colorbar;
            clim([min_max(1), min_max(2)]);
            title([ROI_list{rr} ' ' hemis{hh} ' tstat map: ' active_contrast_list{1}]);

            ROI_table.act_pass_tstat_diff = ROI_table.tstat_act - ROI_table.tstat_pass;
            figure; scatter3(ROI_table, 'R', 'A', 'S', 'filled', 'ColorVariable', 'act_pass_tstat_diff'); colorbar;
            title([ROI_list{rr} ' ' hemis{hh} ' tstat map: (' active_contrast_list{1} ') - (' passive_contrast_list_plot{1} ')']);

            ROI_table = ROI_groupavg_tstats_MD{rr-10,hh,2};
            min_max = [min([ROI_table.tstat_pass; ROI_table.tstat_act]), max([ROI_table.tstat_pass; ROI_table.tstat_act])];

            figure; 
            scatter3(ROI_table,'R', 'A', 'S', 'filled', 'ColorVariable', 'tstat_pass'); colorbar;
            clim([min_max(1), min_max(2)]);
            title([ROI_list{rr} ' ' hemis{hh} ' tstat map: ' passive_contrast_list_plot{4}]);

            figure; scatter3(ROI_table, 'R', 'A', 'S', 'filled', 'ColorVariable', 'tstat_act'); colorbar;
            clim([min_max(1), min_max(2)]);
            title([ROI_list{rr} ' ' hemis{hh} ' tstat map: ' active_contrast_list{4}]);

            ROI_table.act_pass_tstat_diff = ROI_table.tstat_act - ROI_table.tstat_pass;
            figure; scatter3(ROI_table, 'R', 'A', 'S', 'filled', 'ColorVariable', 'act_pass_tstat_diff'); colorbar;
            title([ROI_list{rr} ' ' hemis{hh} ' tstat map: (' active_contrast_list{4} ') - (' passive_contrast_list_plot{4} ')']);

            keyboard;
        end
    end
end


