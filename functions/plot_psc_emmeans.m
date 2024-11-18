function h = plot_psc_emmeans(emmeans_table, title_txt)
%PLOT_PSC_EMMEANS 

%% TODO: make xticklabels the task and xlabel the ROI type in else statement

if nargin == 1
    title_txt = '';
end

set(groot, 'DefaultAxesTickLabelInterpreter', 'none')

h = figure;
if height(emmeans_table) ~= 6
    xs = 1:height(emmeans_table);
    xticklab_inds = 1:height(emmeans_table);
    xtick = 1:height(emmeans_table);
    xlims = [0, height(emmeans_table)+1];
    xticklabs = {emmeans_table.Row{xticklab_inds}};
    xlab = 'Conditions';
else
    % make sure row order of table is as expected
    emmeans_table = sortrows(emmeans_table,'RowNames');
    xs = [1,2,4,5,7,8];
    xticklab_inds = [1,2,3,4,5,6];
    xlims = [0,9];
    xtick = [1,2,4,5,7,8]; 
    xticklabs = cellfun(@(x) [x ' task'], cellstr(emmeans_table{xticklab_inds,2}), 'UniformOutput',false);
    xlab = 'MD ROIs                         Auditory ROIs                         Visual ROIs';
end

% errorbar(xs([1,2,3,4]), emmeans_table.Estimated_Marginal_Mean([1,2,3,4]), emmeans_table.SE([1,2,3,4]), '.', 'MarkerSize', 40, "LineStyle", "none", ...
%     "Color", [0 0.4470 0.7410]); hold on;
% errorbar(xs([5,6,7,8]), emmeans_table.Estimated_Marginal_Mean([5,6,7,8]), emmeans_table.SE([5,6,7,8]), '.', 'MarkerSize', 40, "LineStyle", "none", ...
%     "Color", [0.8500 0.3250 0.0980]);
% errorbar(xs([9,10,11,12]), emmeans_table.Estimated_Marginal_Mean([9,10,11,12]), emmeans_table.SE([9,10,11,12]), '.', 'MarkerSize', 40, "LineStyle", "none", ...
%     "Color", [0.4660 0.6740 0.1880]);
errorbar(xs, emmeans_table.Estimated_Marginal_Mean, emmeans_table.SE, '.', 'MarkerSize', 40, "LineStyle", "none")
hold on;
xlabel(xlab);
xticks(xtick);
xticklabels(xticklabs);
xtickangle(45);
ylabel('Estimated Marginal Mean PSC');
xlim(xlims);
title(title_txt);
grid on;

end

