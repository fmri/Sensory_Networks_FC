%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to convert the BIDs formatted.tsv condition
% timing files to freesurfer formatted .para condition timing files
%
% Tom Possidente July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccc;

%% Set up subject info data
experiment_name = 'spacetime';
task = 'spacetime'; % localizer or spacetime
tsv_fname = 'f_events_taskswitch_PST.tsv'; % f_events.tsv
para_fname = 'taskswitch_spacetime_conditions_PST.para'; % [task '_condition_timing.para']
projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

if strcmp(task, 'localizer')
    seq_name = 'x1WayLocalizer';
    dataDir = [projectDir '/data/unpacked_data_nii_fs_localizer/'];
    fsd = 'localizer';
elseif strcmp(task, 'spacetime')
    seq_name = 'spacetime';
    dataDir = [projectDir '/data/unpacked_data_nii/'];
    fsd = 'bold';
else
    error(['task variable not recognized. Should be "localizer" or "spacetime", instead it is ' task]);
end

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'RR'}));
N = length(subjCodes);

% Tom Localizer condition numbers:
% 1 = Fixation
% 2 = Passive Auditory
% 3 = Passive Tactile
% 4 = Passive Visual
% 5 = Active Auditory
% 6 = Active Tactile
% 7 = Active Visual

% DB localizer
% 1 = vA
% 2 = vP
% 3 = aA
% 4 = aP
% 5 = tA
% 6 = tP
% 7 = f

toms_order = [1,2,3,4,5,6,7];
db_order = [7,4,6,2,3,5,1];

%% Loop over subjs

for ss = 1:N

    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    
    % Get run numbers
    runs = subjDf_cut.([seq_name 'Runs']){subjRow};

    if isempty(runs) && strcmp(task, 'localizer')
        runs = subjDf_cut.('x3WayLocalizerRuns'){subjRow};
    end

    if contains(runs, '/')
        runs = replace(runs, '/', ',');
    end
    runs = str2num(runs);
    
    % Loop over runs
    for rr = 1:length(runs)
        runDir = [dataDir subjCode '/' fsd '/00' num2str(rr) '/'];

        % Load .tsv file
        tsv = readtable([runDir tsv_fname], "FileType","text");

        % Convert to freesurfer para format
        mat = tsv{:,:}; % get rid of header row 
        para = [mat(:,1), mat(:,3), mat(:,2)]; % switch columns 2 and 3
        
        if strcmp(task, 'localizer')
            para(:,2) = db_order(para(:,2));
        end

        % Save as .para
        writematrix(para, [runDir para_fname], 'filetype','text', 'delimiter','\t')

    end

end