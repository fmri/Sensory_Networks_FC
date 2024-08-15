%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to get and plot information on which subjects are
% missing which ROIs and how many vertices are in each ROI
%
% Tom Possidente August 9 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'))
ccc;

%% Get subjIDs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'AH', 'SL'})); % no ROIs for these subjs
N = length(subjCodes);

%% Set key paths and variables
ROI_path = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';

fsaverage = true; % false for individual subj space ROI info

all_ROIs = {'sPCS', 'iPCS', 'midIFS', 'aIPS', 'pIPS', 'DO', 'LOT', 'VOT', ...
    'tgPCS', 'cIFSG', 'CO', 'FO', 'pAud', 'cmSFG', ...
    'aINS', 'preSMA', 'dACC', 'ppreCun'};
ROI_names_lh = cellfun(@(x) [x ' lh'], all_ROIs, 'UniformOutput',false);
ROI_names_rh = cellfun(@(x) [x ' rh'], all_ROIs, 'UniformOutput',false);
ROI_names_lhrh = [ROI_names_lh, ROI_names_rh];
N_ROIs = length(all_ROIs);

%% Loop through subjs and get all ROI information
all_ROI_files = {dir(ROI_path).name};
if fsaverage
    all_ROI_files = all_ROI_files(~contains(all_ROI_files, 'indiv'));
else
    all_ROI_files = all_ROI_files(contains(all_ROI_files, 'indiv'));
end

ROI_matrix = zeros(N,N_ROIs*2);
ROI_size_matrix = nan(N,N_ROIs*2);

for ss = 1:N

    subjID = subjCodes{ss};
    label_fnames = all_ROI_files(contains(all_ROI_files, [subjID '_']));
    if strcmp(subjID, 'NS')
        non_aINS = ~contains(label_fnames, 'aINS');
        good_aINS = contains(label_fnames, 'NS_aINS');
        label_fnames = label_fnames(non_aINS | good_aINS);
    end

    %% left hemisphere
    label_fnames_lh = label_fnames(contains(label_fnames, '_lh'));
    label_names_lh = cellfun(@(x) split(x,'_'), label_fnames_lh, 'UniformOutput', false); % split file name to isolate ROI name
    label_names_lh = string(cellfun(@(x) x{find(contains(x,subjID))+1}, label_names_lh, 'UniformOutput', false)); % ROI name should be 2nd to last item
    ROI_matrix(ss, 1:N_ROIs) = ismember(all_ROIs, label_names_lh);

    for rr = 1:length(label_fnames_lh)
        ROI_ind = find(ismember(all_ROIs, label_names_lh{rr}));
        ROI_size_matrix(ss,ROI_ind) = height(readtable([ROI_path label_fnames_lh{rr}], 'FileType','text'));
    end

    %% right hemisphere
    label_fnames_rh = label_fnames(contains(label_fnames, '_rh'));
    label_names_rh = cellfun(@(x) split(x,'_'), label_fnames_rh, 'UniformOutput', false); % split file name to isolate ROI name
    label_names_rh = string(cellfun(@(x) x{find(contains(x,subjID))+1}, label_names_rh, 'UniformOutput', false)); % ROI name should be 2nd to last item
    ROI_matrix(ss, N_ROIs+1:end) = ismember(all_ROIs, label_names_rh);

    for rr = 1:length(label_fnames_rh)
        ROI_ind = find(ismember(all_ROIs, label_names_rh{rr}));
        ROI_size_matrix(ss,ROI_ind+N_ROIs) = height(readtable([ROI_path label_fnames_rh{rr}], 'FileType','text'));
    end

end

ROI_table = array2table(ROI_matrix, 'VariableNames', ROI_names_lhrh', 'RowNames', subjCodes);
ROI_size_table = array2table(ROI_size_matrix, 'VariableNames', ROI_names_lhrh', 'RowNames', subjCodes);
perc_missing_subj = array2table(sum(ROI_matrix,2) ./ (N_ROIs*2), 'RowNames', subjCodes);
perc_missing_ROI = array2table(sum(ROI_matrix,1)' ./ N, 'RowNames', ROI_names_lhrh);

ROI_matrix_combine_lhrh = [ROI_matrix(:,1:N_ROIs); ROI_matrix(:,(N_ROIs+1):end)];
perc_missing_ROI_lhrh = array2table(sum(ROI_matrix_combine_lhrh,1)' ./ (N*2), 'RowNames', all_ROIs);

%% Make variables necessary for plotting
ROI_matrix_combine_lhrh = [ROI_matrix(:,1:N_ROIs); ROI_matrix(:,(N_ROIs+1):end)];
perc_missing_ROI_lhrh = array2table(sum(ROI_matrix_combine_lhrh,1)' ./ (N*2), 'RowNames', all_ROIs);

ROI_size_mat_lhrh_order_inds_pre = [1:N_ROIs; N_ROIs+1:N_ROIs*2];
ROI_size_mat_lhrh_order_inds = ROI_size_mat_lhrh_order_inds_pre(:)';
ROI_size_mat_lhrh_order = ROI_size_matrix(:,ROI_size_mat_lhrh_order_inds);
ROI_names_lhrh_order = ROI_names_lhrh(ROI_size_mat_lhrh_order_inds);
per_missing_ROI_lhrh_order = perc_missing_ROI.(1)(ROI_size_mat_lhrh_order_inds);


%% Plotting

% Plot overaly ROI size distribution
figure;
histogram(ROI_size_matrix(:), 500)
xlabel('Number of Vertices');
ylabel('Number of ROIs');
title('All ROIs Sizes');

figure;
histogram(ROI_size_matrix(:), 500)
set(gca, 'xscale', 'log');
xlabel('Number of Vertices');
ylabel('Number of ROIs');
title('All ROIs Sizes');

% Scatter plot of ROI size by ROI
figure;
ROI_size_matrix_combine_lhrh = [ROI_size_matrix(:,1:N_ROIs); ROI_size_matrix(:,(N_ROIs+1):end)];
ROI_size_matrix_combine_lhrh = ROI_size_matrix_combine_lhrh(:);
x = repmat(1:N_ROIs, [N*2,1]);
x = x(:);
c = [repmat([1:N]',[1,N_ROIs]); repmat([1:N]',[1,N_ROIs])];
c = c(:);
plt = swarmchart(x,ROI_size_matrix_combine_lhrh, 50, c, 'filled', '');
colormap jet;
xticks(1:N_ROIs);
new_xticklabs = arrayfun(@(x) [all_ROIs{x} ' (' num2str(round(perc_missing_ROI_lhrh.(1)(x)*100, 1)) '% present)'], 1:length(all_ROIs), 'UniformOutput', false);
xticklabels(new_xticklabs); % make tick labels show % missing per ROI also
xlabel('ROI (lh and rh)');
ylabel('Number of Vertices');

% Plot ROI size hemispheric lateralization
figure;
y = ROI_size_mat_lhrh_order(:);
x = repmat(1:N_ROIs*2, [N,1]);
x = x(:);
c = repmat([1:N]',[1,N_ROIs*2]);
c = c(:);
plt = swarmchart(x,y, 50, c, 'filled', '');
colormap jet;
xticks(1:N_ROIs*2);
new_xticklabs = arrayfun(@(x) [ROI_names_lhrh_order{x} ' (' num2str(round(per_missing_ROI_lhrh_order(x)*100, 1)) '%)'], 1:N_ROIs*2, 'UniformOutput', false);
xticklabels(new_xticklabs); % make tick labels show % missing per ROI also
xlabel('ROI (lh and rh)');
ylabel('Number of Vertices');







