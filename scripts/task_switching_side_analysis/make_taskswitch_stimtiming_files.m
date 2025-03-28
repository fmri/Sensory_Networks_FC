%%%%
% The purpose of this script is to create new tsv files for condition
% timing that correspond to task switching times
%
% Created: Tom Possidente - January 2025
%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Get subj codes
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.spacetimeRuns,''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'RR'}));
N = length(subjCodes);

%% Set key variables
csv_dir = '/projectnb/somerslab/vaibhavt/Projects/Trifloc/Data/SubjectsBehavioralFiles/';
ordered_conditions = {'Block Fixation', 'Auditory Passive', 'Tactile Passive', 'Visual Passive',...
    'Auditory Spatial', 'Tactile Spatial', 'Visual Spatial', 'Auditory Temporal',...
    'Tactile Temporal', 'Visual Temporal'}; 
condition_names = {'fixation', 'switch_passive', 'switch_spatial', 'switch_temporal'};
N_conditions = length(condition_names);
tsv_pathbase = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii/';

% 1 = fixation
% 2 = aP = auditory passive
% 3 = tP = tactile passive
% 4 = vP = visual passive
% 5 = aS = auditory spatial
% 6 = tS = tactile spatial
% 7 = vS = visual spatial
% 8 = aT = auditory temporal
% 9 = tT = tactile temporal
% 10 = vT = visual temporal

%% Extract condition info from csvs for each run for each subj

for ss = 1:N

    % Get file names of stim timing csvs
    subjCode = subjCodes{ss};
    row_ind = ismember(subjDf_cut.subjCode, subjCode);
    spacetime_date = subjDf_cut.spacetimeDate{row_ind};
    N_runs = length(split(subjDf_cut.spacetimeRuns{row_ind}, [",", "/"]));
    subj_stim_files = {dir([csv_dir spacetime_date subjCode]).name};
    file_mask = (contains(subj_stim_files, {'spatial_temporal'}) | contains(subj_stim_files, {'stt_ac_'})) & contains(subj_stim_files, {'.csv'}) & ~contains(subj_stim_files, '.~lock');
    subj_stim_files = subj_stim_files(file_mask);
    assert(length(subj_stim_files)==N_runs, ['Subj ' subjCode ': found ' num2str(length(subj_stim_files)) ' behavioral files but subj has ' num2str(N_runs) ' runs']);

    [subj_stim_files, order_correct] = check_stimfile_order(subj_stim_files);
    assert(order_correct, 'Stim files are not in sequential time order from earliest to latest filename timestamp');

    % Loop over runs and extract condition timing information
    for rr = 1:N_runs

        beh_data = readtable([csv_dir spacetime_date subjCode '/' subj_stim_files{rr}], 'Delimiter', ',');
        modalities = beh_data.modality;
        types = beh_data.type;
        trial_onsets = beh_data.fixation_block_started;
        cue_started = beh_data.cue_text_started;
        N_rows = height(beh_data);

        conditions = nan(N_rows-1,1);
        for row = 1:N_rows-1 % last row always empty
            condition_str = [modalities{row} ' ' types{row}];
            conditions(row) = find(ismember(ordered_conditions, condition_str));
        end

        % Chunk conditions to remove consecutive repeats
        [conditions, n_per_condition, change_detect] = chunk(conditions, true);
        assert(all(conditions<11), 'One or more condition label is unrecognized (greater than 10)');
        assert(length(unique(conditions))==10, 'All 10 conditions are not represented');
        
        % Get associated onset times
        onsets = trial_onsets(change_detect);
        onsets(isnan(onsets)) = cue_started(find(strcmpi(types,'Fixation'),1,'first')) + 3; % fixation does not have a trial_onset value, so use cue_started and add 3s
        onsets = onsets - (onsets(1)-3); % for condition files, time starts at first cue time (3s before first trial fixation onset)

        % Check to make sure all passive/fixation have 3 consecutive occurences (30s) and the rest have 6 (60s)
        % and make duration column at the same time
        durations = nan(length(conditions),1);
        for cc = 1:length(conditions)
            if ismember(conditions(cc), [1,2,3,4])
                assert(ismember(n_per_condition(cc), [3,4]) , 'Unexpected length of condition');
                durations(cc) = 30;
            else
                assert(ismember(n_per_condition(cc), [6,8]), 'Unexpected length of condition');
                durations(cc) = 60;
            end
        end
        
        % Add interblock conditions for task switching analysis 
        [onsets_ts, durations_ts, conditions_ts] = make_taskswitch_conditions(onsets, durations, conditions, ordered_conditions, condition_names, 1:length(condition_names));
        
        % Actually make the tsv files
        condition_table_taskswitch = array2table([onsets_ts, durations_ts, conditions_ts], 'VariableNames', {'onset', 'duration', 'trial_type'});
        writetable(condition_table_taskswitch, [tsv_pathbase, subjCode, '/bold/00', num2str(rr) '/f_events_taskswitch_PST.tsv'], 'Delimiter', '\t', 'FileType','text');
    end

    disp(['Finished subj ' subjCode])

end



%% Helper function(s)

function [stim_files, order_correct] = check_stimfile_order(stim_files)
% Makes sure the order of the stimulus files is sequential in time by
% checking the timestamps present in the file names. These need to be in
% order so that they correctly correspond to runs 1,2,3,4... 

    suffixes = cellfun(@(x) x(end-15:end-4), stim_files, 'UniformOutput',false);
    
    if contains(suffixes{1}, '.') % there are 2 different possible formats for the timestamps in the filenames
        suffixes = cellfun(@(x) split(x(end-11:end-7),'h'), suffixes, 'UniformOutput',false);
        times = cellfun(@(x) str2double(x{1})+(str2double(x{2})/60), suffixes, 'UniformOutput',false); % converting to hours    
    elseif contains(suffixes{1}, '_')
        suffixes = cellfun(@(x) split(x(end-6:end),'_'), suffixes, 'UniformOutput',false);
        day = cellfun(@(x) x{1}, suffixes, 'UniformOutput',false);
        same_day = all(strcmp(day, day{1})); % all should be the same day
        assert(same_day, 'Not all stim files have the same day in their file name timestamp');
        times = cellfun(@(x) str2double(x{2}), suffixes, 'UniformOutput',false);
    else
        error(['Unexpected file name format for subj ' subjCode ': ' suffixes{1}]);% 1 = fixation
% 2 = aP = auditory passive
% 3 = tP = tactile passive
% 4 = vP = visual passive
% 5 = aS = auditory spatial
% 6 = tS = tactile spatial
% 7 = vS = visual spatial
% 8 = aT = auditory temporal
% 9 = tT = tactile temporal
% 10 = vT = visual temporal
    end
    
    for tt = 1:length(times)-1
        assert(times{tt} < times{tt+1}, 'Stim files are not in sequential time order from earliest to latest filename timestamp');
    end
    order_correct = true;

end

function [onsets_out, durations_out, conditions_out] = make_taskswitch_conditions(onsets, durations, conditions, ordered_conditions, taskswitch_conditions, taskswitch_condition_labels)
% 
    
    %% Identify which task switch is happening between each condition and get its label number
    conditions_out = [];
    for ii = 1:length(conditions)
        curr_cond_name = ordered_conditions{conditions(ii)};
        if ismember(curr_cond_name, {'Block Fixation'})
            conditions_out = [conditions_out; 1];
            continue;
        elseif ismember(curr_cond_name, {'Auditory Passive', 'Visual Passive', 'Tactile Passive'})
            curr_cond_type = 'switch_passive';
        elseif ismember(curr_cond_name, {'Visual Spatial', 'Auditory Spatial', 'Tactile Spatial'})
            curr_cond_type = 'switch_spatial';
        elseif ismember(curr_cond_name, {'Visual Temporal',  'Tactile Temporal', 'Auditory Temporal'})
            curr_cond_type = 'switch_temporal';
        else
            error('condition name not recognized');
        end
        conditions_out = [conditions_out; taskswitch_condition_labels(ismember(taskswitch_conditions, curr_cond_type))];
    end
    durations_out = repelem(3,length(conditions_out),1); % all task switch durations are 3s
    durations_out(conditions_out==1) = durations(conditions_out==1);
    onsets_out = onsets - 3;

end






