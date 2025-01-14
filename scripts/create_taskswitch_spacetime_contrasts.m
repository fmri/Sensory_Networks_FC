%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to make task switching contrasts out of the
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
subjCodes = subjCodes(~ismember(subjCodes, {'RR'}));
N = length(subjCodes);

fsaverage = true;

data_dir = [projectDir 'data/unpacked_data_nii/'];
analysis_name_lh = 'taskswitch_VAT_spacetime_contrasts_lh_polyfit2hrf1';
analysis_name_rh = 'taskswitch_VAT_spacetime_contrasts_rh_polyfit2hrf1';
TR = 2;
design_type = 'blocked';
para_name = 'taskswitch_spacetime_conditions_VAT.para';
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

refeventdur = '3';
nconditions = '4'; %%%%%%%%%%%%

% spacetime condition orders
% 1 = fixation
% 2 = switch to visual block
% 3 = switch to auditory block
% 4 = switch to tactile block

%% Create analysis/contrast info

% run make analysis lh/rh
% using -event-related because freesurfer uses this to mean either event-related or block design. Using polyfit 1 (can use 2 if
% we want more noise regressors). spmhrf 0 to use 0 derivatives of the hrf.
% mcextreg option will use 12 motion regressors (3 translation 3 rotation
% and derivatives)

if run_mkanalysis

    unix(['mkanalysis-sess -a ' analysis_name_lh ' -funcstem ' funcstem_lh ...
        ' -surface ' space ' lh -fsd bold -event-related -paradigm ' para_name ...
        ' -nconditions ' nconditions ' -refeventdur ' refeventdur ' -TR ' num2str(TR) ...
        ' -polyfit 2 -spmhrf 1 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])

    unix(['mkanalysis-sess -a ' analysis_name_rh ' -funcstem ' funcstem_rh ...
        ' -surface ' space ' rh -fsd bold -event-related -paradigm ' para_name ...
        ' -nconditions ' nconditions ' -refeventdur ' refeventdur ' -TR ' num2str(TR) ...
        ' -polyfit 2 -spmhrf 1 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])
end


% run make contrast
if run_mkcontrast
    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'ts_visual-f -a 2 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'ts_visual-f -a 2 -c 1'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'ts_auditory-f -a 3 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'ts_auditory-f -a 3 -c 1'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'ts_tactile-f -a 4 -c 1'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'ts_tactile-f -a 4 -c 1'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'ts_visual-ts_audtact -a 2 -c 3 -c 4'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'ts_visual-ts_audtact -a 2 -c 3 -c 4'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'ts_auditory-ts_vistact -a 3 -c 2 -c 4'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'ts_auditory-ts_vistact -a 3 -c 2 -c 4'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'ts_tactile-ts_visaud -a 4 -c 2 -c 3'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'ts_tactile-ts_visaud -a 4 -c 2 -c 3'])
end



%% Loop through subjs
for ss = 1:N

    subjCode = subjCodes{ss};

    % run glm
    unix(['selxavg3-sess -s ' subjCode ' -d ' data_dir ' -analysis ' analysis_name_lh...
        ' -no-preproc -overwrite'])
    unix(['selxavg3-sess -s ' subjCode ' -d ' data_dir ' -analysis ' analysis_name_rh...
        ' -no-preproc -overwrite'])

end


