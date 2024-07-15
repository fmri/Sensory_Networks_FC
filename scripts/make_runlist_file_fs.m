%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to
%%% Tom Possidente - July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
N = height(subjDf_cut);
rlf_name = 'localizer_contrasts_runlistfile.txt';

for ss = 1:N
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));

    % Get functional scan run numbers
    seq_name = 'x1WayLocalizer';
    func_runs = subjDf_cut.([seq_name, 'Runs']){subjRow};
    if isempty(func_runs)
        seq_name = 'x3WayLocalizer';
        func_runs = subjDf_cut.([seq_name, 'Runs']){subjRow};
    end

    if contains(func_runs, '/') % runs with different fieldmaps
        func_runs = replace(func_runs, '/', ','); % still take all runs
    end
    func_runs = str2num(func_runs);

    % create run list file 
    nruns = length(func_runs);
    unix(['for i in $(seq 1 ' num2str(nruns) '); do echo "00$i"; done >> ' ...
        projectDir '/data/unpacked_data_nii_fs_localizer/' subjCode '/localizer/' rlf_name])

end




