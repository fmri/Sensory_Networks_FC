%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to make task contrasts out of the
%%% freesurfer preprocessed spacetime task data surfaces.
%%% Tom Possidente - September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'RR', 'AH', 'SL'}));
N = length(subjCodes);

fsaverage = true;

data_dir = [projectDir 'data/unpacked_data_nii/'];
analysis_name_lh = 'spacetime_contrasts_lh_xtra_contrasts_NoyceRep';
analysis_name_rh = 'spacetime_contrasts_rh_xtra_contrasts_NoyceRep';
TR = 2;
design_type = 'blocked';
para_name = 'spacetime_condition_timing_tom.para';
if fsaverage
    space = 'fsaverage';
else
    space = 'self';
end
funcstem_lh = ['fmcpr_tu.siemens.sm3.' space '.lh.nii.gz'];
funcstem_rh = ['fmcpr_tu.siemens.sm3.' space '.rh.nii.gz'];
rlf_name = 'spacetime_contrasts_runlistfile.txt';

run_mkanalysis = true;
run_mkcontrast = true;

refeventdur = '30';
nconditions = '10';

% spacetime condition orders
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




%% Create analysis/contrast info

% run make analysis lh/rh
% using -event-related even though this is a block design because this option only matters for
% preprocessing which has already been done. Using polyfit 1 (can use 2 if
% we want more noise regressors). spmhrf 0 to use 0 derivatives of the hrf.
% mcextreg option will use 12 motion regressors (3 translation 3 rotation
% and derivatives)

if run_mkanalysis
    unix(['mkanalysis-sess -a ' analysis_name_lh ' -funcstem ' funcstem_lh ...
        ' -surface ' space ' lh -fsd bold -event-related -paradigm ' para_name ...
        ' -nconditions ' nconditions ' -refeventdur ' refeventdur ' -TR ' num2str(TR) ...
        ' -polyfit 1 -spmhrf 0 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])

    unix(['mkanalysis-sess -a ' analysis_name_rh ' -funcstem ' funcstem_rh ...
        ' -surface ' space ' rh -fsd bold -event-related -paradigm ' para_name ...
        ' -nconditions ' nconditions ' -refeventdur ' refeventdur ' -TR ' num2str(TR) ...
        ' -polyfit 1 -spmhrf 0 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])
end


% run make contrast
if run_mkcontrast

    % visual passive - auditory passive
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'vP-aP -a 4 -c 2'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'vP-aP -a 4 -c 2'])

    % auditory passive - visual passive
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'aP-vP -a 2 -c 4'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'aP-vP -a 2 -c 4'])

end



%% Loop through subjs
parfor ss = 1:N

    subjCode = subjCodes{ss};

    % run glm
    unix(['selxavg3-sess -s ' subjCode ' -d ' data_dir ' -analysis ' analysis_name_lh...
        ' -no-preproc -overwrite'])
    unix(['selxavg3-sess -s ' subjCode ' -d ' data_dir ' -analysis ' analysis_name_rh...
        ' -no-preproc -overwrite'])

end


