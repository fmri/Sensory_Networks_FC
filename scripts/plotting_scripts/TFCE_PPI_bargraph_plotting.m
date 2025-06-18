%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to plot the effect sizes from the TFCE
% gPPI analyses comparing connectivity changes between auditory and visual
% WM tasks 
% Tom Possidente - March 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/'));
ccc;

%% Load effect size data
load('/projectnb/somerslab/tom/projects/sensory_networks_FC/data/conn_toolbox_folder/conn_localizer_task/results/secondlevel/avsm_gPPI/AllSubjects/aA(1).vA(-1)/effect_sizes.mat', 'data', 'data_maxCI', 'data_minCI', 'names');
es_data_avsm = data([2,1],:)';
es_max_avsm = data_maxCI([2,1],:)';
es_min_avsm = data_minCI([2,1],:)';

%% Plot
clusters_avsm = {'cmSFG(R), iPCS-sm(R), midIFS(R)', 'midIFS(R), iPCS-sm(R), pVis', 'FO(R), midIFS, iPCS-sm(R)',...
                'cIFSG(L), preSMA-sm(L), aINS-sm(L)', 'cIFSG, FO, cmSFG(L)', 'FO(R), preSMA-sm(R), aINS-sm(R)',...
                'pAud, iPCS-sm(L), midIFS/MFG-sm(L)', 'iPCS-sm(L), pVis', 'sPCS-sm, sPCS',...
                'cIFSG(L), sPCS', 'iPCS(L), sPCS(L)', 'cIFSG(L), sPCS-sm(R)'};
figure; 
hb = bar(es_data_avsm); % get the bar handles
hold on;
for k = 1:size(es_data_avsm,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, es_data_avsm(:,k), es_min_avsm(:,k)-es_data_avsm(:,k), es_max_avsm(:,k)-es_data_avsm(:,k), 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1);
end

% Set Axis properties
set(gca,'xticklabel', clusters_avsm);
ylabel('Effect Size')
xlabel('ROI Clusters');
legend({'Visual WM', 'Auditory WM'});
title('Change in connectivity during visual and auditory working memory')
set(gca, 'FontSize', 18)
