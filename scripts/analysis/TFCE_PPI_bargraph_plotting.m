%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to plot the effect sizes from the TFCE
% gPPI analyses comparing connectivity changes between auditory and visual
% WM tasks 
% Tom Possidente - March 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Load effect size data
% load('/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_localizer_task/results/secondlevel/gPPI_3sm_21/AllSubjects/aA(1).vA(-1)/effect_size.mat', 'data', 'data_maxCI', 'data_minCI', 'names');
% es_data = data([2,1],:)';
% es_max = data_maxCI([2,1],:)';
% es_min = data_minCI([2,1],:)';
% load('/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_localizer_task/results/secondlevel/gPPI_3sm_21/AllSubjects/aA(1).vA(-1)/effect_size_supramodal.mat', 'data', 'data_maxCI', 'data_minCI', 'names');
% es_data_sm = data([2,1])';
% es_max_sm = data_maxCI([2,1])';
% es_min_sm = data_minCI([2,1])';
load('/projectnb/somerslab/tom/projects/spacetime_network/data/conn_toolbox_folder/conn_localizer_task/results/secondlevel/avsm_gPPI/AllSubjects/aA(1).vA(-1)/effect_sizes.mat', 'data', 'data_maxCI', 'data_minCI', 'names');
es_data_avsm = data([2,1],:)';
es_max_avsm = data_maxCI([2,1],:)';
es_min_avsm = data_minCI([2,1],:)';

%% Plot
clusters_avsm = {'cmSFG(R), sm_iPCS(R), midIFS(R)', 'midIFS(R), sm_iPCS(R), pVis', 'FO(R), midIFS, sm_iPCS(R)',...
                'cIFSG(L), sm_preSMA(L), sm_aINS(L)', 'cIFSG, FO, cmSFG(L)', 'FO(R), sm_preSMA(R), sm_aINS(R)',...
                'pAud, sm_iPCS(L), sm_midFSG(L)', 'sm_iPCS(L), pVis', 'sm_sPCS, sPCS',...
                'cIFSG(L), sPCS', 'iPCS(L), sPCS(L)', 'cIFSG(L), sm_sPCS(R)'};
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

% %% Plot
% clusters_sm = {'sPCS, midIFS, pVis', 'cIFS/G, CO (L), cmSFG (L), tgPCS (R), pAud', 'cIFS/G, cmSFG, FO', 'FO, preSMA', 'preSMA, dACC (L), pAud',...
%             'iPCS (R), pVis', 'tgPCS (R), pAud', 'cIFS/G, preSMA', 'tgPCS (L), aINS, dACC (R)'};
% figure; 
% hb = bar(es_data); % get the bar handles
% hold on;
% for k = 1:size(es_data,2)
%     % get x positions per group
%     xpos = hb(k).XData + hb(k).XOffset;
%     % draw errorbar
%     errorbar(xpos, es_data(:,k), es_min(:,k)-es_data(:,k), es_max(:,k)-es_data(:,k), 'LineStyle', 'none', ...
%         'Color', 'k', 'LineWidth', 1);
% end
% 
% % Set Axis properties
% set(gca,'xticklabel', clusters_sm);
% ylim([-0.2 0.4])
% ylabel('Effect Size')
% xlabel('ROI Clusters');
% grid on;
% legend({'Visual WM', 'Auditory WM'});
% title('Change in connectivity during visual and auditory working memory')
% set(gca, 'FontSize', 18)
% 
% %% Plot
% clusters_sm = {'sm_iPCS, pVis'};
% figure; 
% 
% hb1 = bar(1, es_data_sm(1)); hold on;
% hb2 = bar(2, es_data_sm(2)); hold on;
% errorbar([1,2], es_data_sm, es_min_sm-es_data_sm, es_max_sm-es_data_sm, 'LineStyle', 'none', ...
%     'Color', 'k', 'LineWidth', 1);
% 
% % Set Axis properties
% set(gca, 'XTick', 1.5)
% set(gca,'xticklabel', clusters_sm, 'TickLabelInterpreter', 'none');
% ylim([-0.25 0.1])
% ylabel('Effect Size')
% xlabel('ROI Clusters');
% grid on;
% legend({'Visual WM', 'Auditory WM'});
% title('Change in connectivity during visual and auditory working memory')
% set(gca, 'FontSize', 18)
