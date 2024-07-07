%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to make sure I am using the same run
% numbers for localizer and spacetime tasks as Vaihbav did.
%
% Tom Possidente July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccc;

%% Set up subject info data
experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
N = length(subjCodes);

vaibhav_dirbase = '/projectnb/somerslab/trifloc2/'; 

if strcmp(experiment_name, 'spacetime')
    vaibhav_filename = 'SpatialTemporal.txt';
elseif strcmp(experiment_name, 'x1WayLocalizer')
    vaibhav_filename = '1wayPilot.txt';
else
    error('Unknown experiment name, valid choices are "spacetime" and "x1WayLocalizer"');
end

%% Loop over subjs

for ss = 1:N

    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    experiment_date = subjDf_cut.([experiment_name,'Date']){subjRow};
    
    % Get Tom run numbers
    tom_runs = subjDf_cut.([experiment_name 'Runs']){subjRow};

    if strcmp(experiment_name, 'x1WayLocalizer') && isempty(tom_runs)
        tom_runs = subjDf_cut.('3WayLocalizerRuns'){subjRow};
    end

    if contains(tom_runs, '/')
        tom_runs = replace(tom_runs, '/', ',');
    end
    tom_runs = str2num(tom_runs);
    
    % Get Vaibhav run numbers
    vaibhav_dir = [vaibhav_dirbase experiment_date subjCode '/' vaibhav_filename];
    if strcmp(experiment_name, 'x1WayLocalizer') && ~isfile(vaibhav_dir)
        vaibhav_dir = [vaibhav_dirbase experiment_date subjCode '/1way3wayPilot.txt'];
    end
    vaibhav_runs_table = readtable(vaibhav_dir);
    vaibhav_runs = vaibhav_runs_table{:,1};

    % Compare and display results
    issame = ismember(tom_runs, vaibhav_runs);

    if all(issame)
        disp(['Subj ' subjCode ' ' experiment_name ': same run numbers'])
        disp(' ');
    else
        disp(['Subj ' subjCode ' ' experiment_name ': DIFFERENT RUN NUMEBRS']);
        disp(['Toms runs: ' num2str(tom_runs)]);
        disp(['Vaibhav runs: ' num2str(vaibhav_runs')]);
        disp(' ');
    end

end


