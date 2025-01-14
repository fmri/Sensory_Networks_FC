%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to make task switching contrasts out of the
%%% freesurfer preprocessed spacetime task data surfaces.
%%% Tom Possidente - January 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Set up directories and subj info

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'RR', 'SL'}));
subjCodes = {'AH'}
N = length(subjCodes);

fsaverage = true;

data_dir = [projectDir 'data/unpacked_data_nii/'];
analysis_name_lh = 'spacetime_contrasts_lh_newcondfiles';
analysis_name_rh = 'spacetime_contrasts_rh_newcondfiles';
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

    % unix(['mkanalysis-sess -a ' analysis_name_lh ' -funcstem ' funcstem_lh ...
    %     ' -surface ' space ' lh -fsd bold -event-related -paradigm ' para_name ...
    %     ' -nconditions ' nconditions ' -refeventdur ' refeventdur ' -TR ' num2str(TR) ...
    %     ' -polyfit 2 -spmhrf 1 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])
    % 
    % unix(['mkanalysis-sess -a ' analysis_name_rh ' -funcstem ' funcstem_rh ...
    %     ' -surface ' space ' rh -fsd bold -event-related -paradigm ' para_name ...
    %     ' -nconditions ' nconditions ' -refeventdur ' refeventdur ' -TR ' num2str(TR) ...
    %     ' -polyfit 2 -spmhrf 1 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])
end


% run make contrast
if run_mkcontrast
    % All active - passive
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'A-P -a 5 -a 7 -a 8 -a 10 -a 6 -a 9 -c 2 -c 4 -c 3'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'A-P -a 5 -a 7 -a 8 -a 10 -a 6 -a 9 -c 2 -c 4 -c 3'])

    % visual and audio active - passive
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'AV_A-P -a 5 -a 7 -a 8 -a 10 -c 2 -c 4'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'AV_A-P -a 5 -a 7 -a 8 -a 10 -c 2 -c 4'])

    % single modality all active - passive
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sVtV-pV -a 7 -a 10 -c 4'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sVtV-pV -a 7 -a 10 -c 4'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sAtA-pA -a 5 -a 8 -c 2'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sAtA-pA -a 5 -a 8 -c 2'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sTtT-pT -a 6 -a 9 -c 3'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sTtT-pT -a 6 -a 9 -c 3'])

    % single modality spatial active - passive
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sV-pV -a 7 -c 4'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sV-pV -a 7 -c 4'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sA-pA -a 5 -c 2'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sA-pA -a 5 -c 2'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sT-pT -a 6 -c 3'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sT-pT -a 6 -c 3'])

    % single modality temporal active - passive
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'tV-pV -a 10 -c 4'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'tV-pV -a 10 -c 4'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'tA-pA -a 8 -c 2'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'tA-pA -a 8 -c 2'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'tT-pT -a 9 -c 3'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'tT-pT -a 9 -c 3'])

    % Active modality - acitve modality
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sVtV-sAtA -a 7 -a 10 -c 5 -c 8'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sVtV-sAtA -a 7 -a 10 -c 5 -c 8'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sVtV-sTtT -a 7 -a 10 -c 6 -c 9'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sVtV-sTtT -a 7 -a 10 -c 6 -c 9'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sAtA-sTtT -a 5 -a 8 -c 6 -c 9'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sVtV-sTtT -a 7 -a 10 -c 6 -c 9'])

    % Passive modality - fixation
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'pV-f -a 4 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'pV-f -a 4 -c 1'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'pA-f -a 2 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'pA-f -a 2 -c 1'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'pT-f -a 3 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'pT-f -a 3 -c 1'])

    % Active spatial modality - fixation
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sV-f -a 7 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sV-f -a 7 -c 1'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sA-f -a 5 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sA-f -a 5 -c 1'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'sT-f -a 6 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'sT-f -a 6 -c 1'])

    % Active temporal modality - fixation
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'tV-f -a 10 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'tV-f -a 10 -c 1'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'tA-f -a 8 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'tA-f -a 8 -c 1'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'tT-f -a 9 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'tT-f -a 9 -c 1'])


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


