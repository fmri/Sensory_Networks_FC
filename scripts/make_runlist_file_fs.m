%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to
%%% Tom Possidente - July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';
task = 'spacetime'; % localizer or spacetime
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
N = height(subjDf_cut);
rlf_name = [task '_contrasts_runlistfile.txt'];

% for ss = 1:N
%     subjCode = subjCodes{ss};
%     subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
% 
%     % Get functional scan run numbers
%     func_runs = subjDf_cut.([seq_name, 'Runs']){subjRow};
%     if isempty(func_runs) && strcmp(task, 'localizer')
%         seq_name_alt = 'x3WayLocalizer';
%         func_runs = subjDf_cut.([seq_name_alt, 'Runs']){subjRow};
%     end
% 
%     if contains(func_runs, '/') % runs with different fieldmaps
%         func_runs = replace(func_runs, '/', ','); % still take all runs
%     end
%     func_runs = str2num(func_runs);
% 
%     % create run list file 
%     nruns = length(func_runs);
%     unix(['for i in $(seq 1 ' num2str(nruns) '); do echo "00$i"; done >> ' ...
%         dataDir subjCode '/' fsd '/' rlf_name])
% 
% end



%% If you want to create separate run list files for each run:
for ss = 1:N
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));

    % Get functional scan run numbers
    func_runs = subjDf_cut.([seq_name, 'Runs']){subjRow};
    if isempty(func_runs) && strcmp(task, 'localizer')
        seq_name_alt = 'x3WayLocalizer';
        func_runs = subjDf_cut.([seq_name_alt, 'Runs']){subjRow};
    end

    if contains(func_runs, '/') % runs with different fieldmaps
        func_runs = replace(func_runs, '/', ','); % still take all runs
    end
    func_runs = str2num(func_runs);

    % create run list file 
    nruns = length(func_runs);
    for rr = 1:nruns
        rlf_name = [task '_contrasts_runlistfile' num2str(rr) '.txt'];
        unix(['echo "00' num2str(rr) '" >> ' dataDir subjCode '/' fsd '/' rlf_name])
    end

end

