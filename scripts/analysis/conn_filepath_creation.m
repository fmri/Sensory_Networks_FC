%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to print to console the filepaths for all
% T1s, fieldmaps, and functional runs, so that they can be copy/pasted into
% Conn Toolbox for easy initialization of data directories
%
% Created: Tom Possidente - March 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('/projectnb/somerslab/tom/helper_functions/');
addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

supramodal_ROIs = false; % give annot paths to supramodal ROI set instead of full vbias, abias, posterior, MD set
func_topupapplied = true; % give paths for functionals with fmaps already applied
resting_state = false;
localizer = true;

% Load in subject info
subjDf = load_subjInfo();
if resting_state 
    data_dir = 'rest';
    experiment_name = 'rest';
    subjInds = ~strcmp(subjDf.([experiment_name,'Runs']),'');
    reject_subjs = {'MM','RR','AH','PQ','RT','SL'}; % MM rs in different dims, PQ no subthreshold motion runs, RR low activation, AH no ROIs due to no localizer runs w/o motion, RT no rs runs at all
elseif localizer
    data_dir = 'localizer';
    experiment_name = 'x1WayLocalizer';
    subjInds = ~( strcmp(subjDf.([experiment_name,'Runs']),'') & strcmp(subjDf.('x3WayLocalizerRuns'),'') );
    reject_subjs = {'RR', 'AH', 'SL'};
else
    data_dir = 'bold';
    experiment_name = 'spacetime';
    subjInds = ~strcmp(subjDf.([experiment_name,'Runs']),'');
    reject_subjs = {'RR', 'AH', 'SL', 'AI'};
end
subjDf_cut = subjDf(subjInds,:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, reject_subjs));

ROI_path = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs/';
struct_path = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/';
func_path = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii/';

% Print ROI annot paths for ROI selection
not_found = zeros(length(subjCodes),1);
for ss=1:length(subjCodes)
    if supramodal_ROIs
        subjROIpath = [ROI_path '/lh.' subjCodes{ss} '_sm_ROIs.annot'];
    else
        subjROIpath = [ROI_path '/lh.' subjCodes{ss} '_ROIs.annot'];
    end
    if ~isfile(subjROIpath)
        not_found(ss) = 1;
    else
        disp(subjROIpath)
    end
end
disp(['No ROIs: ' string(subjCodes(logical(not_found)))' ])

% print ROI path
for ss=1:length(subjCodes)
    if supramodal_ROIs
        subjROIpath = [ROI_path subjCodes{ss} '_smROIs.surf.nii'];
    else
        subjROIpath = [ROI_path subjCodes{ss} '_avsm_ROIs.surf.nii'];
    end
    disp(subjROIpath)
end


% Print T1 paths
for ss=1:length(subjCodes)

    subjDirStruct = [struct_path '/' subjCodes{ss} '/mri/T1.nii'];
    assert(isfile(subjDirStruct), ['Subj ' subjCodes{ss} ' T1 file not found'])
    disp(subjDirStruct)

end

% Print functional paths
runs_all = {};
if func_topupapplied
    prefix = 'u';
    suffix = '_topupApplied.nii';
else
    prefix = '';
    suffix = '.nii';
end
for ss=1:length(subjCodes)
    subjCode = subjCodes{ss};
    subjRow = find(strcmp(subjDf_cut.subjCode, subjCode));
    runs = subjDf_cut.([experiment_name, 'Runs']){subjRow};

    if localizer && isempty(runs)
        runs = subjDf_cut.('x3WayLocalizerRuns'){subjRow};
    end

    if contains(runs, '/') % runs with different fieldmaps
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);
    runs_all{ss} = runs;
    for ii=1:length(runs)
        subjDirFunc = [func_path '/' subjCode '/' data_dir, '/00' num2str(ii) '/' prefix 'f' suffix];
        if ~isfile(subjDirFunc)
            subjDirFunc = [subjDirFunc '.gz'];
            if ~isfile(subjDirFunc)
                error(['Subj ' subjCode ' run ' num2str(runs(ii)) ' functional file not found'])
            end
        end
        disp(subjDirFunc)
    end

end

% Print realignment files (rp_f.txt)
for ss=1:length(subjCodes)
    subjCode = subjCodes{ss};
    for ii=1:length(runs_all{ss})
        realignment_filepath = [func_path '/' subjCode '/' data_dir '/00' num2str(ii) '/rp_f.txt'];
        assert(isfile(realignment_filepath), ['Subj ' subjCode ' run ' num2str(ii) ' realignment file not found'])
        disp(realignment_filepath)
    end
end

% Print number of runs per subj
run_nums = nan(length(subjCodes), 1);
for ss = 1:length(subjCodes)
    disp(num2str(length(runs_all{ss})));
    run_nums(ss) = length(runs_all{ss});
end

run_nums_tbl = table(subjCodes, run_nums);
