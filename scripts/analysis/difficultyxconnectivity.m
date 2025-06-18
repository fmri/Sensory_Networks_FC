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
lm_table_vis.Properties.VariableNames = {'subject', 'beta_diff' 'vis_difficulty'};
lm_table_aud.Properties.VariableNames = {'subject', 'beta_diff', 'aud_difficulty'};

lm_table_aud_plot = lm_table_aud;
lm_table_vis_plot = lm_table_vis;
lm_table_vis.beta_diff = lm_table_vis.beta_diff - mean(lm_table_vis.beta_diff); 
lm_table_aud.beta_diff = lm_table_aud.beta_diff - mean(lm_table_aud.beta_diff); 

lm_aud = fitlm(lm_table_aud, 'aud_difficulty ~ 1 + beta_diff')
lm_aud_plot = fitlm(lm_table_aud_plot, 'aud_difficulty ~ 1 + beta_diff');
f = figure; p = plot(lm_aud_plot); title(''); xlabel('Auditory Connections: Auditory Betas - Visual Betas'); ylabel('Auditory WM Change in Frequency');

lm_vis = fitlm(lm_table_vis, 'vis_difficulty ~ 1 + beta_diff')
lm_vis_plot = fitlm(lm_table_vis_plot, 'vis_difficulty ~ 1 + beta_diff');
figure; plot(lm_vis_plot); title(''); xlabel('Visual Connections: Auditory Betas - Visual Betas'); ylabel('Visual WM Change in Frequency');

%% Test individual connections
for ss = 1:N
    % Get network-level changes in connectivity for subj
    subj_data = data_table(data_table.subject==ss,:);
    for cc = 1:N_conns
        % Vis betas
        betas(ss,cc,1) = mean(subj_data.vis_beta(ismember(subj_data.connection_type, connections{cc})));

        % Aud betas
        betas(ss,cc,2) = mean(subj_data.aud_beta(ismember(subj_data.connection_type, connections{cc})));
    end
end

%% Plot
for cc = 1:N_conns
    % vis_tbl = table(betas(:,cc,1), avg_difficulty(:,3));
    % vis_tbl.Properties.VariableNames = {'betas', 'vis_difficulty'};
    % vis_tbl.vis_difficulty = 60 - vis_tbl.vis_difficulty;
    % lm_vis = fitlm(vis_tbl, 'vis_difficulty ~ 1 + betas');
    % figure; plot(lm_vis);
    % xlabel('Visual Betas'); ylabel('Visual WM Task Performance'); 
    % title(['Visual WM: ' connections{cc} ' | p=' num2str(lm_vis.Coefficients{2,4})]);
    % set(gca, 'FontSize', 18);
    % 
    % vis_tbl = table(betas(:,cc,2), avg_difficulty(:,3));
    % vis_tbl.Properties.VariableNames = {'betas', 'vis_difficulty'};
    % vis_tbl.vis_difficulty = 60 - vis_tbl.vis_difficulty;
    % lm_vis = fitlm(vis_tbl, 'vis_difficulty ~ 1 + betas');
    % figure; plot(lm_vis);
    % xlabel('Auditory Betas'); ylabel('Visual WM Task Performance'); 
    % title(['Auditory WM: ' connections{cc} ' | p=' num2str(lm_vis.Coefficients{2,4})]);
    
    aud_tbl = table(betas(:,cc,1), avg_difficulty(:,1));
    aud_tbl.Properties.VariableNames = {'betas', 'aud_difficulty'};
    lm_aud = fitlm(aud_tbl, 'aud_difficulty ~ 1 + betas')
    figure; plot(lm_aud);
    xlabel('Visual Betas'); ylabel('Auditory WM Task Performance'); 
    title(['Visual WM: ' connections{cc} ' | p=' num2str(lm_aud.Coefficients{2,4})]);
    set(gca, 'FontSize', 18);

    aud_tbl = table(betas(:,cc,2), avg_difficulty(:,1));
    aud_tbl.Properties.VariableNames = {'betas', 'aud_difficulty'};
    lm_aud = fitlm(aud_tbl, 'aud_difficulty ~ 1 + betas')
    figure; plot(lm_aud);
    xlabel('Auditory Betas'); ylabel('Auditory WM Task Performance'); 
    title(['Auditory WM: ' connections{cc} ' | p=' num2str(lm_aud.Coefficients{2,4})]);
    set(gca, 'FontSize', 18);

end

% tbl = array2table(betas(:,:,2));
% tbl.Properties.VariableNames = strrep(connections, '<->', '');
% tbl.aud_difficulty = avg_difficulty(:,1);
% fitlm(tbl, 'aud_difficulty ~ 1 + abiasabias + abiasvbias + abiassupramodal + pVissupramodal + pVisvbias')
