%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to retrieve and organize trifloc/spacetime
% localizer stimulus timing data from the timing files created by
% connectome workbench to tsv files that can be read in easily by Conn
% Toolbox
%
% Created: Tom Possidente - July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Set up directories and key variables
projectDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/';
stimTimingBase = '/projectnb/somerslab/hcp_pipeline_subjects/';
experiment_name = 'x1WayLocalizer';
targetDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/behavioral/stim_timing/';

use_3way = false; 

if use_3way
    file_suffix = '3WayLocalizer';
    EV_fileNames = {'3way_Active_Auditory.txt', '3way_Active_Tactile.txt', '3way_Active_Visual.txt',...
        '3way_Fixation_Block.txt', '3way_Passive_Visual.txt'};
else
    file_suffix = '1WayLocalizer';
    EV_fileNames = {'1way_Fixation_Block.txt', '1way_Passive_Auditory.txt', '1way_Passive_Tactile.txt', ...
    '1way_Passive_Visual.txt', '1way_Active_Auditory.txt', '1way_Active_Tactile.txt'...
    '1way_Active_Visual.txt'};
end

condition_codes = 1:length(EV_fileNames);



%% Load subject info
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),'') | ~strcmp(subjDf.('x3WayLocalizerRuns'),''),:);
subjCodes = subjDf_cut.subjCode;
n = length(subjCodes);

%% Loop through subjects
for ss = 1:n

    subjCode = subjCodes{ss};
    runDataDir = [stimTimingBase, lower(subjCode), '3p20/MNINonLinear/Results/'];
    folders = {dir(runDataDir).name};

    % Get number of localizer runs (-1 for folder with metadata)
    filestr = '1wayPilot';
    num_localizer_runs = sum(contains(folders, filestr)) - 1;

    if num_localizer_runs<=0 || (strcmp(subjCode, 'MM') && use_3way) % some subjs have it labeled differently
        filestr = '1way3wayPilot';
        num_localizer_runs = sum(contains(folders, filestr)) - 1;
    end

    % Loop through spacetime runs and extract timing for each
    count = 0;
    for rr = 1:num_localizer_runs

        timing_path = [runDataDir, filestr, num2str(rr), '/EVs/'];
        timing_data = table();
        if use_3way
            num_stimtiming_files = sum(contains({dir(timing_path).name}, '3way'));
            condition = num_stimtiming_files == 5;     
        else
            num_stimtiming_files = sum(contains({dir(timing_path).name}, '1way'));
            condition = num_stimtiming_files == 6 || num_stimtiming_files == 7; % should contain 6 or 7 files (some don't have fixation block)
        end

        if condition
            count = count + 1;
            if ~isfile([targetDir, subjCode, '_run', num2str(count) '_' file_suffix '.txt'])
                for tt = 1:length(EV_fileNames)
                    if isfile([timing_path, EV_fileNames{tt}])
                        time_tbl = readtable([timing_path, EV_fileNames{tt}]);
                        time_tbl(:,3) = num2cell(repelem(condition_codes(tt), height(time_tbl))');
                        timing_data = [timing_data; time_tbl];
                    end
                end

                timing_data.Properties.VariableNames = ["onset","duration","trial_type"];
                timing_data = sortrows(timing_data, 'onset');

                % Save to tsv file
                writetable(timing_data, [targetDir, subjCode, '_run', num2str(count) '_' file_suffix], 'Delimiter', '\t', 'FileType','text');

            end
        end
    end

end
