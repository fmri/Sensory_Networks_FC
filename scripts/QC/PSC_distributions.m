%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to quality check the distributions of PSCs for each
% ROI across all subjs for each condition to look for outliers
% Tom Possidente - November 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load PSC data
load('/projectnb/somerslab/tom/projects/spacetime_network/data/PSCs/PSC_results.mat', 'ROIs','contrasts','psc_results','subjCodes');

%% Plot distributions
%for rr = 1:length(ROIs)
    for cc = 1:length(contrasts)-2
        data = psc_results(:,:,cc,:);
        figure;
        histogram(data(:),20);
        %title([ROIs{rr}, ' ' contrasts{cc}]);
        title([contrasts{cc}, ' mean=' num2str(nanmean(data(:))), ' sd=', num2str(nanstd(data(:)))]);
    end
%end