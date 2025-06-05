%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to extract the staircasing difficulty
%%% levels for each trifloc WM task for each subj and average them over
%%% runs to get an idea of task aptitude/performance for each subj/run
%%% Tom Possidente - May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Initialize Key Variables
%reject_subjs = {'AH', 'SL', 'RR'};
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N = length(subjCodes);

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.('spacetimeRuns'),''),:);
taskname = 'x1WayLocalizer';

behavior_pathbase = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/behavioral/behavioral/';
modalities = {'Auditory', 'Tactile', 'Visual'};

difficulty = cell(N,1);

%% Loop through subjs, extract behavioral files, parse files, get staircased difficulty
for ss = 1:N
    subjCode = subjCodes{ss};
    row_ind = strcmp(subjCode,subjDf_cut.subjCode);
    date = subjDf_cut.([taskname 'Date']){row_ind};
    runs = subjDf_cut.([taskname, 'Runs']){row_ind};
    if contains(runs, '/') 
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);
    num_runs = length(runs);
    subjID = [date upper(subjCode)];

    files = {dir([behavior_pathbase subjID]).name};
    if isempty(files)
        error(['No behavrioal data files for subj ' subjCode]);
    end
    
    files_csv = files( contains(files, '.csv') & (contains(files, 'trifloc_task_') | contains(files, 'ac_trifloc_task')) ); % Get only csv files in this dir

    assert(length(files_csv)==num_runs, ['Subj ' subjCode ': Number of runs specified in subjInfo.csv does not match number of runs found in behavioral data files']);
    
    %difficulty_subj = nan(3, num_runs, 16);
    difficulty_subj = nan(3, num_runs, 16);

    % Loop through each run and collect difficulty data for each task
    for rr = 1:num_runs
        behavioral_data = readtable([behavior_pathbase subjID, '/', files_csv{rr}], 'Delimiter',',');
        for ii = 1:3
            mask = ismember(behavioral_data.modality, modalities{ii}) & ismember(behavioral_data.type, 'Active');
            task_data = behavioral_data(mask,:);
            %difficulty_subj(ii,rr,:) = task_data.trials_intensity;
            difficulty_subj(ii,rr,:,1) = task_data.odd_trial_changeValue;
            difficulty_subj(ii,rr,:,2) = task_data.even_trial_changeValue;
        end
    end
    difficulty{ss} = difficulty_subj;
end

%% Compute averages
avg_difficulty = nan(N,3);

for ss = 1:N
    for ii = 1:3
        %avg_difficulty(ss,ii) = mean(difficulty{ss}(ii,:,:), [2,3])
        avg_difficulty(ss,ii) = squeeze(mean(difficulty{ss}(ii,:,:,:), [2,3,4]));
    end
end

%% Plot averages
figure;
histogram(avg_difficulty(:,1), 'NumBins',5);
hold;
histogram(avg_difficulty(:,2), 'NumBins',5);
histogram(avg_difficulty(:,3), 'NumBins',5);

save('staircase_difficulties.mat', 'difficulty', 'avg_difficulty', 'subjCodes', 'modalities');


