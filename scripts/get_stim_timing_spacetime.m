%%%%
% The purpose of this script is to retrieve and organize trifloc/spacetime
% task stimulus timing data from the timing files created by
% connectome workbench to tsv files that can be read in easily by Conn
% Toolbox
%
% Created: Tom Possidente - Feb 2024
%%%%%

addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Set up directories and key variables
projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
stimTimingBase = '/projectnb/somerslab/hcp_pipeline_subjects/';
experiment_name = 'spacetime';
EV_fileNames = {'Fixation_Block.txt', 'Passive_Auditory.txt', 'Passive_Tactile.txt', 'Passive_Visual.txt',...
    'Spatial_Auditory.txt', 'Spatial_Tactile.txt', 'Spatial_Visual.txt',...
    'Temporal_Auditory.txt', 'Temporal_Tactile.txt', 'Temporal_Visual.txt'};
condition_codes = 1:length(EV_fileNames);
targetDir = '/projectnb/somerslab/tom/projects/spacetime_network/data/behavioral/stim_timing/';


%% Load subject info
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = {'AF'}
n = length(subjCodes);

%% Loop through subjects
for ss = 1:n

    subjCode = subjCodes{ss};
    runDataDir = [stimTimingBase, lower(subjCode), '3p20/MNINonLinear/Results/'];
    folders = {dir(runDataDir).name};

    % Get number of spatial temporal runs (-1 for folder with metadata)
    num_spacetime_runs = sum(contains(folders, 'SpatialTemporal')) - 1;

    % Loop through spacetime runs and extract timing for each
    for rr = 1:num_spacetime_runs

        %if ~isfile([targetDir, subjCode, '_run', num2str(rr) '.txt'])
            timing_path = [runDataDir, 'SpatialTemporal', num2str(rr), '/EVs/'];
            timing_data = nan(length(condition_codes), 3);
            for tt = 1:length(EV_fileNames)
                time_tbl = readtable([timing_path, EV_fileNames{tt}]);
                timing_data(tt,:) = time_tbl{1,:};
                timing_data(tt,3) = condition_codes(tt);
            end

            timing_datatbl = array2table(timing_data, 'VariableNames', ["onset","duration","trial_type"]);
            timing_datatbl = sortrows(timing_datatbl, 'onset');

            % Save to tsv file
            writetable(timing_datatbl, [targetDir, subjCode, '_run', num2str(rr)], 'Delimiter', '\t', 'FileType','text');
        %end

    end

end
