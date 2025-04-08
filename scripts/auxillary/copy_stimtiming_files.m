%%%%
% The purpose of this script is to copy the relevant stimulus timing files
% to each subject's functional scan directory so that Conn Toolbox can find
% it automatically.
%
% Created: Tom Possidente - March 2024
%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir = [projectDir 'data/'];
stimtimeDir = [subjectsDir 'behavioral/stim_timing/'];

experiment_name = 'x1WayLocalizer'; % Change stimtiming files to x1WayLocalizer here if you want to copy those files into the functional localizer dirs

if strcmp(experiment_name, 'spacetime')
    file_suffix = '.txt';
    func_dir = 'bold';
elseif strcmp(experiment_name, 'x1WayLocalizer')
    file_suffix = '_1Waylocalizer.txt';
    func_dir = 'localizer';
elseif strcmp(experiment_name, 'x3WayLocalizer')
    file_suffix = '_3WayLocalizer.txt';
    func_dir = 'localizer';
else
    error('Experiment name not recognized. Should be either "spacetime" or "x1WayLocalizer"');
end

%%

for ss = 1:length(subjCodes)

    %% Set up subj varibles
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));

    % Get run numbers of func data
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};
    if isempty(runs)
        runs = subjDf_cut.('x3WayLocalizerRuns'){subjRow};
    end
    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);

    for rr = 1:length(runs)
        sourceFile = [stimtimeDir subjCode '_run' num2str(rr) file_suffix];
        destFile = [subjectsDir 'unpacked_data_nii_fs_localizer/' subjCode '/' func_dir '/00' num2str(rr) '/f_events.tsv'];
        unix(['cp ' sourceFile ' ' destFile])
    end


end








