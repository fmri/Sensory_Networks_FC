%%%%
% The purpose of this script is to load in the behavioral data for the
% spacetime localizer task and calculate the % correct for each subj for each
% condition for each modality and run statistical testing on them
%
% Created: Tom Possidente - Jan 2025
%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

save_out = false;

%% Load subject info
taskname = 'x1WayLocalizer';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.('spacetimeRuns'),''),:);
subjCodes = subjDf_cut.subjCode;
n = length(subjCodes);

%% Loop through subjs and get % correct
behavioral_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/behavioral/behavioral/';
modality_names = {'visual', 'auditory', 'tactile'};
n_conditions = 2*length(modality_names);
correct_colnames = {'trials_response', 'trialsresponse', 'odd_trial_responseCorrect'};
perc_correct_all = nan(n,n_conditions);
n_cond = nan(n,n_conditions);
perc_correct_byrun = nan(n,2,6);

if strcmp(taskname, 'x1WayLocalizer')
    condition_str = {'active', 'passive'};
    condition_order = {'visual_active', 'visual_passive', 'auditory_active', 'auditory_passive', 'tactile_active', 'tactile_passive'};
elseif strcmp(taskname, 'spacetime')
    condition_str = {'spatial', 'temporal'};
    condition_order = {'visual_spatial', 'visual_temporal', 'auditory_spatial', 'auditory_temporal', 'tactile_spatial', 'tactile_temporal'};
end

for ss = 1:n

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

    files = {dir([behavioral_dir subjID]).name};
    if isempty(files)
        disp(['No behavrioal data files for subj ' subjCode]);
        continue
    end

    if strcmp(taskname, 'spacetime')
        files_csv = files( contains(files, '.csv') & (contains(files, '_spatial_temporal_') | contains(files, 'ac_trifloc_task')) ); % Get only csv files in this dir
    elseif strcmp(taskname, 'x1WayLocalizer')
        files_csv = files( contains(files, '.csv') & (contains(files, 'trifloc_task_') | contains(files, 'ac_trifloc_task')) ); % Get only csv files in this dir
    else
        error('unrecognized taskname');
    end

    assert(length(files_csv)==num_runs, ['Subj ' subjCode ': Number of runs specified in subjInfo.csv does not match number of runs found in behavioral data files']);

    responses = [];
    modalities = {};
    conditions = {};


    % Loop through each run and collect reponses, modality, and condition data
    len_runs = zeros(num_runs+1,1);
    for rr = 1:num_runs
        behavioral_data = readtable([behavioral_dir subjID, '/', files_csv{rr}], 'Delimiter',',');
        nrows = height(behavioral_data);
        modalities(ss,rr,:) = behavioral_data.modality;
        conditions(ss,rr,:) = behavioral_data.type;

        which_column = find(ismember(correct_colnames, behavioral_data.Properties.VariableNames), 1);
        if which_column ~= 1
            keyboard;
        end
        responses(ss,rr,:) = behavioral_data.(correct_colnames{which_column});
        len_runs(rr+1) = length(behavioral_data.(correct_colnames{which_column}));
    end

    subj_responses = squeeze(responses(ss,:,:));

    % Loop through each modality and calculate % correct
    for mm = 1:length(modality_names)
        modality_cond1_mask = squeeze(strcmpi(modality_names{mm}, modalities(ss,:,:)) & strcmpi(condition_str{1}, conditions(ss,:,:)));
        n_cond(ss,(mm*2)-1) = sum(modality_cond1_mask,'all');
        perc_correct = mean(subj_responses(modality_cond1_mask), 'all', 'omitnan');
        perc_correct_all(ss,(mm*2)-1) = perc_correct;
        assert(~isnan(perc_correct_all(ss,(mm*2)-1)), 'percent correct calculated as nan, should not happen');

        if perc_correct <= 0.55
            disp(['Subj ' subjCode ' condition ' modality_names{mm} ' cond1 has behavior performance below 55% (' num2str(perc_correct) ')']);
        end

        for rr = 1:num_runs
            perc_correct_byrun(ss,mm,rr) = mean(responses(ss,rr, modality_cond1_mask(rr,:)), 'omitnan');
        end

        modality_cond2_mask = squeeze(strcmpi(modality_names{mm}, modalities(ss,:,:)) & strcmpi(condition_str{2}, conditions(ss,:,:)));
        n_cond(ss,(mm*2)) = sum(modality_cond2_mask, 'all');
        perc_correct = mean(subj_responses(modality_cond2_mask), 'all', 'omitnan');
        perc_correct_all(ss,(mm*2)) = perc_correct;
        assert(~isnan(perc_correct_all(ss,(mm*2))) || strcmp(condition_str{2}, 'passive'), 'percent correct calculated as nan, should not happen');
        
        if perc_correct <= 0.55
            disp(['Subj ' subjCode ' condition ' modality_names{mm} ' cond2 has behavior performance below 55% (' num2str(perc_correct) ')']);
        end

    end

end

%% Save out results
perc_correct_all = array2table(perc_correct_all(:,[1,3,5]), 'VariableNames',condition_order([1,3,5]))
if save_out
    save('behavioral_percent_correct_data.mat', 'subjCodes', 'perc_correct_all', 'n_cond');
    writetable(perc_correct_all, 'behavioral_data.csv')
end

missing_all_data = all(ismissing(perc_correct_all),2);
perc_correct_all=perc_correct_all(~missing_all_data,:); % delete any rows with all nans
n = n-sum(missing_all_data);

%% Power analysis
% For a repeated measures ANOVA With sample N=20, alpha = 0.05, power = 0.8,
% and a 3x2 measurements, we need to calculate the correlation among
% repeated measures and assess nonsphericity in order to do a power
% analysis.

% Calculating correlation among repeated measures
perc_correct_matrix = table2array(perc_correct_all);
corr_matrix = triu(corr(perc_correct_matrix),1); % get upper triangle (without diag) of correlation matrix
mean_corr = mean(corr_matrix(corr_matrix~=0),'all', 'omitnan'); % average correlation
% This process is not technically statistically valid, but it gives a rough
% idea of the correlation among measures, good enough for a power analysis
% % We should have fisher z transformed, averaged, convert back to
% correlation coeff.

% Looking at sphericity: assumption of equal variances in the differences
% between measures
measure_differences = nan(n, 6); % with 5 measures, there are 10 comparisons
count = 0;
for ii = 1:3
    for jj = 1:3
        if (ii ~= jj) && (ii > jj)
            count = count + 1;
            measure_differences(:,count) = perc_correct_matrix(:,ii) - perc_correct_matrix(:,jj);
            disp(num2str([ii, jj])); % checking all combos are done
        end

    end
end

vars = var(measure_differences); % if these variances are similar, sphericity is good

% gpower was used to calculate effect size detectable from a reapeated
% measures, within factors ANOVA with the following parameters:
% alpha: 0.05
% power = 0.8
% total sample size = 24
% number of groups = 1
% number of measurements = 3
% corr among repeated measures = 0.4286
% nonspereicity correction (epsilon) = 0.8
% The result is an effect size detectable of f = 0.30997 which is
% approximately a cohen's d of 0.620 (medium-large effect size)

% To translate a cohens d of 0.620 into a raw mean difference we can use the
% fact that cohens d = difference in means / pooled SD
pooled_SD = sqrt(mean(sqrt(var(perc_correct_matrix))));
diff_in_means_detectable = 0.620*pooled_SD;

% We can detect a difference in means of about 0.1864 as a main effect
% currently. Interaction effects likely require a larger difference in
% means

%%
LME_design = perc_correct_all(:,[1,3,5]);
%LME_design = LME_design(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}),:);
LME_design.subject = [1:n]';


%% Actually run the 2 way repeated measures ANOVA
design_tbl = table([0,1,2]', 'VariableNames', {'modality'});
design_tbl.modality = categorical(design_tbl.modality);
rm = fitrm(perc_correct_all, "visual_active,auditory_active,tactile_active~1", WithinDesign=design_tbl);
results = ranova(rm,'WithinModel','modality')


%% Plot swarmplot for visualization
close all;
figure;
boxplot(perc_correct_matrix, 'Labels', condition_order([1,3,5]));
ylim([0,1]);
ylabel('Proportion Correct');

newmat_modality = [mean(perc_correct_matrix(:,[1,2]),2), mean(perc_correct_matrix(:,[3,4]),2), mean(perc_correct_matrix(:,[5,6]),2)];
newmat_task = [mean(perc_correct_matrix(:,[1,3,5]),2), mean(perc_correct_matrix(:,[2,4,6]),2)];

figure;
boxplot(newmat_modality, 'Labels', {'visual', 'auditory', 'tactile'});
ylabel('Proportion Correct');
ylim([0,1]);

figure;
boxplot(newmat_task, 'Labels', {'spatial', 'temporal'});
ylabel('Proportion Correct');
ylim([0,1]);



