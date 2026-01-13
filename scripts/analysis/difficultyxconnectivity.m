%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to investigate correlations between
%%% behavioral performance and connectivity between networks for each WM
%%% modality 
%%% Tom Possidente - May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Initialize Key Variables
%reject_subjs = {'AH', 'SL', 'RR'};
%subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};

% Load difficulties
load('/projectnb/somerslab/tom/projects/sensory_networks_FC/data/behavioral/staircase_difficulties_oddevenavg.mat', 'avg_difficulty', 'modalities', 'subjCodes');
diff_subjCodes = subjCodes;

% Load gppi betas 
load('/projectnb/somerslab/tom/projects/sensory_networks_FC/data/gPPI_va-aa_fulltable.mat', 'data_table', 'subjCodes')
assert(all(strcmp(diff_subjCodes, subjCodes)))
N = length(subjCodes);

connections = {'pVis<->supramodal', 'pVis<->vbias', 'abias<->abias', 'abias<->vbias', 'abias<->supramodal'};
N_conns = length(connections);

betas = nan(N,N_conns,2);

%% Group into aud and vis hypotheses
connection_types = {{'abias<->abias', 'abias<->vbias', 'abias<->supramodal'}, {'pVis<->supramodal', 'pVis<->vbias', }};
betas_grouped = nan(N, length(connection_types));
lm_table_vis = table();
lm_table_aud = table();

for ss = 1:N
    % Get network-level changes in connectivity for subj
    subj_data = data_table(data_table.subject==ss,:);
    for cc = 1:length(connection_types)
        conntype_inds = ismember(subj_data.connection_type, connection_types{cc});
        beta_diffs = subj_data.aud_beta(conntype_inds) - subj_data.vis_beta(conntype_inds);
        betas_grouped(ss,cc) = mean(beta_diffs);

        if cc==1
            lm_table_aud(end+1,:) = {ss, betas_grouped(ss,cc), avg_difficulty(ss,1)};
        else
            lm_table_vis(end+1,:) = {ss, betas_grouped(ss,cc), avg_difficulty(ss,3)};
        end
    end
end

lm_table_vis.Properties.VariableNames = {'subject', 'beta_diff_visconns' 'vis_difficulty'};
lm_table_aud.Properties.VariableNames = {'subject', 'beta_diff_audconns', 'aud_difficulty'};

lm_table_all = table(lm_table_vis.subject, lm_table_vis.beta_diff_visconns, lm_table_aud.beta_diff_audconns, lm_table_vis.vis_difficulty, lm_table_aud.aud_difficulty);
lm_table_all.Properties.VariableNames = {'subject', 'beta_diff_visconns', 'beta_diff_audconns', 'vis_difficulty', 'aud_difficulty'};

lm_aud = fitlm(lm_table_aud, 'aud_difficulty ~ 1 + beta_diff_audconns')
f = figure; p = plot(lm_aud); title(''); xlabel('Auditory Connections: Auditory Betas - Visual Betas'); ylabel('Auditory WM Change in Frequency');

lm_vis = fitlm(lm_table_vis, 'vis_difficulty ~ 1 + beta_diff_visconns')
figure; plot(lm_vis); title(''); xlabel('Visual Connections: Auditory Betas - Visual Betas'); ylabel('Visual WM Change in Frequency');

lm_audbeta_predvis = fitlm(lm_table_all, 'vis_difficulty ~ 1 + beta_diff_audconns')
figure; plot(lm_audbeta_predvis); title(''); xlabel('Auditory Connections: Auditory Betas - Visual Betas'); ylabel('Visual WM Change in Frequency');

lm_visbeta_predaud = fitlm(lm_table_all, 'aud_difficulty ~ 1 + beta_diff_visconns')
figure; plot(lm_visbeta_predaud); title(''); xlabel('Visual Connections: Auditory Betas - Visual Betas'); ylabel('Auditory WM Change in Frequency');