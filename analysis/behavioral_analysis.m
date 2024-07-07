%%%%
% The purpose of this script is to load in the behavioral data for the 
% spacetime task and calculate the % correct for each subj for each
% condition for each modality and run statistical testing on them
%
% Created: Tom Possidente - Feb 2024
%%%%

addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load subject info 
taskname = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([taskname,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
n = length(subjCodes);

%% Loop through subjs and get % correct
behavioral_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/behavioral/behavioral/';
modality_names = {'visual', 'auditory', 'tactile'};
task_names = {'spatial', 'temporal'};
condition_order = {'visual_spatial', 'visual_temporal', 'auditory_spatial', 'auditory_temporal', 'tactile_spatial', 'tactile_temporal'};
n_conditions = length(task_names)*length(modality_names);
perc_correct_all = nan(n,n_conditions);
n_cond = nan(n,n_conditions); 

for ss = 1:n

    subjCode = subjCodes{ss};
    row_ind = strcmp(subjCode,subjCodes);
    spactime_date = subjDf_cut.([taskname 'Date']){row_ind};
    runs = subjDf_cut.([taskname, 'Runs']){row_ind};
    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);
    num_runs = length(runs);
    subjID = [spactime_date upper(subjCode)];

    files = {dir([behavioral_dir subjID]).name};
    if isempty(files)
        disp(['No behavrioal data files for subj ' subjCode]);
        continue
    end
    files_csv = files(contains(files, '.csv')); % Get only csv files in this dir
    assert(length(files_csv)==num_runs, ['Subj ' subjCode ': Number of runs specified in subjInfo.csv does not match number of runs found in behavioral data files']);
    
    responses = [];
    modalities = {};
    conditions = {};

    % Loop through each run and collect reponses, modality, and condition data
    for rr = 1:num_runs 
        behavioral_data = readtable([behavioral_dir subjID, '/', files_csv{rr}], 'Delimiter',',');
        nrows = height(behavioral_data);
        modalities(length(modalities)+1:length(modalities)+nrows) = behavioral_data.modality;
        conditions(length(conditions)+1:length(conditions)+nrows) = behavioral_data.type;
        try
            responses(length(responses)+1:length(responses)+nrows) = behavioral_data.trials_response;
        catch
            disp(['"trials_response" column not found for subj ' subjCode ' run ' num2str(rr) ' using "odd_trial_responseCorrect instead.'])
            responses(length(responses)+1:length(responses)+nrows) = behavioral_data.odd_trial_responseCorrect;
        end
    end

    % Loop through each modality and calculate % correct
    for mm = 1:length(modality_names)
        modality_spatial_mask = strcmpi(modality_names{mm}, modalities) & strcmpi('spatial', conditions);
        n_cond(ss,(mm*2)-1) = sum(modality_spatial_mask);
        perc_correct_all(ss,(mm*2)-1) = mean(responses(modality_spatial_mask));
        assert(~isnan(perc_correct_all(ss,(mm*2)-1)), 'percent correct calculated as nan, should not happen');

        modality_temporal_mask = strcmpi(modality_names{mm}, modalities) & strcmpi('temporal', conditions);
        n_cond(ss,(mm*2)) = sum(modality_temporal_mask);
        perc_correct_all(ss,(mm*2)) = mean(responses(modality_temporal_mask));
        assert(~isnan(perc_correct_all(ss,(mm*2))), 'percent correct calculated as nan, should not happen');
    end

    
end

%% Save out results
perc_correct_all = array2table(perc_correct_all, 'VariableNames',condition_order)
save('behavioral_percent_correct_data.mat', 'subjCodes', 'perc_correct_all', 'n_cond');
writetable(perc_correct_all, 'behavioral_data.csv')

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
mean_corr = mean(corr_matrix(corr_matrix~=0),'all'); % average correlation
% This process is not technically statistically valid, but it gives a rough
% idea of the correlation among measures, good enough for a power analysis
% % We should have fisher z transformed, averaged, convert back to
% correlation coeff.

% Looking at sphericity: assumption of equal variances in the differences
% between measures
measure_differences = nan(n, 10); % with 5 measures, there are 10 comparisons
count = 0;
for ii = 1:5
    for jj = 1:5
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
% number of measurements = 6
% corr among repeated measures = 0.3184
% nonspereicity correction (epsilon) = 0.8 
% The result is an effect size detectable of f = 0.273 which is
% approximately a cohen's d of 0.546 (medium-large effect size)

% To translate a cohens d of 0.546 into a raw mean difference we can use the
% fact that cohens d = difference in means / pooled SD
pooled_SD = sqrt(mean(sqrt(var(perc_correct_matrix))));
diff_in_means_detectable = 0.546*pooled_SD;

% We can detect a difference in means of about 0.1900 as a main effect
% currently. Interaction effects likely require a larger difference in
% means


%% Actually run the 2 way repeated measures ANOVA
design_tbl = table([1,0,1,0,1,0]', [0,0,1,1,2,2]', 'VariableNames', {'task', 'modality'});
design_tbl.task = categorical(design_tbl.task);
design_tbl.modality = categorical(design_tbl.modality);
rm = fitrm(perc_correct_all, "visual_spatial,visual_temporal,auditory_spatial,auditory_temporal,tactile_spatial,tactile_temporal~1", WithinDesign=design_tbl);
results = ranova(rm,'WithinModel','task*modality')

design_tbl = table([1,0,1,0]', [0,0,1,1]', 'VariableNames', {'task', 'modality'});
design_tbl.task = categorical(design_tbl.task);
design_tbl.modality = categorical(design_tbl.modality);
rm = fitrm(perc_correct_all([1:9,11,13:20],1:4), "visual_spatial,visual_temporal,auditory_spatial,auditory_temporal~1", WithinDesign=design_tbl);
results = ranova(rm,'WithinModel','task*modality')

% design_tbl = table([1,0,1,0]', [0,0,1,1]', 'VariableNames', {'task', 'modality'});
% design_tbl.task = categorical(design_tbl.task);
% design_tbl.modality = categorical(design_tbl.modality);
% rm = fitrm(perc_correct_all(:,[1,2,5,6]), "visual_spatial,visual_temporal,tactile_spatial,tactile_temporal~1", WithinDesign=design_tbl);
% results = ranova(rm,'WithinModel','task*modality')
% 
% design_tbl = table([1,0,1,0]', [0,0,1,1]', 'VariableNames', {'task', 'modality'});
% design_tbl.task = categorical(design_tbl.task);
% design_tbl.modality = categorical(design_tbl.modality);
% rm = fitrm(perc_correct_all(:,[3,4,5,6]), "auditory_spatial,auditory_temporal,tactile_spatial,tactile_temporal~1", WithinDesign=design_tbl);
% results = ranova(rm,'WithinModel','task*modality')

%% Plot swarmplot for visualization
close all;
figure; 
boxplot(perc_correct_matrix, 'Labels', condition_order);
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



