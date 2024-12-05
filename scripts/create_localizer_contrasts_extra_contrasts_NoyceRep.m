%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to make task contrasts out of the
%%% freesurfer preprocessed localizer data surfaces - in order to draw ROIs
%%% from the contrast maps
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
N = length(subjCodes);

fsaverage = true;

data_dir = [projectDir 'data/unpacked_data_nii_fs_localizer/'];
analysis_name_lh = 'localizer_contrasts_lh_xtra_contrasts_NoyceRep'; %%%%%%%%%%%%%%
analysis_name_rh = 'localizer_contrasts_rh_xtra_contrasts_NoyceRep'; %%%%%%%%%%%%%%
TR = 2;
design_type = 'blocked';
para_name = 'localizer_condition_timing.para';
if fsaverage
    space = 'fsaverage';
else
    space = 'self';
end
funcstem_lh = ['fmcpr_tu.siemens.sm5.' space '.lh.nii.gz'];
funcstem_rh = ['fmcpr_tu.siemens.sm5.' space '.rh.nii.gz'];
rlf_name = 'localizer_contrasts_runlistfile.txt'; %%%%%%%%%%%%%%%%

run_mkanalysis = true;
run_mkcontrast = true;

use_3way_contrasts = false;

if use_3way_contrasts == true
    nconditions = 5;
else
    nconditions = 7;
end


% localizer condition orders (1way)
% 1 = vA
% 2 = vP
% 3 = aA
% 4 = aP
% 5 = tA
% 6 = tP
% 7 = f

% localizer condition orders (3way)
% 1 = aA
% 2 = tA
% 3 = vA
% 4 = f
% 5 = p


%% Create analysis/contrast info

% run make analysis lh/rh
% using -event-related even
% though this is a block design because this option only matters for
% preprocessing which has already been done. Using polyfit 1 (can use 2 if
% we want more noise regressors). spmhrf 0 to use 0 derivatives of the hrf.
% mcextreg option will use 12 motion regressors (3 translation 3 rotation
% and derivatives)

if run_mkanalysis
    unix(['mkanalysis-sess -a ' analysis_name_lh ' -funcstem ' funcstem_lh ...
        ' -surface ' space ' lh -fsd localizer -event-related -paradigm ' para_name ...
        ' -nconditions ' num2str(nconditions) ' -refeventdur 32 -TR ' num2str(TR) ...
        ' -polyfit 1 -spmhrf 0 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])

    unix(['mkanalysis-sess -a ' analysis_name_rh ' -funcstem ' funcstem_rh ...
        ' -surface ' space ' rh -fsd localizer -event-related -paradigm ' para_name ...
        ' -nconditions ' num2str(nconditions) ' -refeventdur 32 -TR ' num2str(TR) ...
        ' -polyfit 1 -spmhrf 0 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])
end


% run make contrast

% localizer condition orders
% 1 = vA
% 2 = vP
% 3 = aA
% 4 = aP
% 5 = tA
% 6 = tP
% 7 = f

% localizer condition orders (3way)
% 1 = aA
% 2 = tA
% 3 = vA
% 4 = f
% 5 = p

if run_mkcontrast

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'vP-aP -a 4 -c 2'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'vP-aP -a 4 -c 2'])

    unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
        'aP-vP -a 2 -c 4'])
    unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
        'aP-vP -a 2 -c 4'])


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

%%%%%%%%%% Run for subjs RR MM PP (no fixation condition) %%%%%%%%%%

%%
% analysis_name_lh_nofix = 'localizer_contrasts_lh_nofix_run1';
% analysis_name_rh_nofix = 'localizer_contrasts_rh_nofix_run1';
% rlf_name = 'localizer_contrasts_runlistfile_run1.txt';
% subjCodes = {'RR', 'MM', 'PP'};
%
% run_mkanalysis = true;
% run_mkcontrast = true;
%
% if run_mkanalysis
%     unix(['mkanalysis-sess -a ' analysis_name_lh_nofix ' -funcstem ' funcstem_lh ...
%         ' -surface fsaverage lh -fsd localizer -event-related -paradigm ' para_name ...
%         ' -nconditions 6 -refeventdur 32 -TR ' num2str(TR) ...
%         ' -polyfit 1 -spmhrf 0 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])
%
%     unix(['mkanalysis-sess -a ' analysis_name_rh_nofix ' -funcstem ' funcstem_rh ...
%         ' -surface fsaverage rh -fsd localizer -event-related -paradigm ' para_name ...
%         ' -nconditions 6 -refeventdur 32 -TR ' num2str(TR) ...
%         ' -polyfit 1 -spmhrf 0 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])
% end
%
%
% if run_mkcontrast
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'A-P -a 1 -a 3 -a 5 -c 2 -c 4 -c 6'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'A-P -a 1 -a 3 -a 5 -c 2 -c 4 -c 6'])
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'vA-vP -a 1 -c 2'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'vA-vP -a 1 -c 2'])
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'aA-aP -a 3 -c 4'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'aA-aP -a 3 -c 4'])
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'tA-tP -a 5 -c 6'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'tA-tP -a 5 -c 6'])
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'V-A -a 1 -a 2 -c 3 -c 4'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'V-A -a 1 -a 2 -c 3 -c 4'])
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'V-T -a 1 -a 2 -c 5 -c 6'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'V-T -a 1 -a 2 -c 5 -c 6'])
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'A-T -a 3 -a 4 -c 5 -c 6'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'A-T -a 3 -a 4 -c 5 -c 6'])
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'vA-aAtA -a 1 -c 3 -c 5'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'vA-aAtA -a 1 -c 3 -c 5'])
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'aA-vAtA -a 3 -c 1 -c 5'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'aA-vAtA -a 3 -c 1 -c 5'])
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'tA-vAaA -a 5 -c 1 -c 3'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'tA-vAaA -a 5 -c 1 -c 3'])
%
%     unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
%         'vAaA-vPaP -a 1 -a 3 -c 2 -c 4'])
%     unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
%         'vAaA-vPaP -a 1 -a 3 -c 2 -c 4'])
%
% end
%
% parfor ss = 1:length(subjCodes)
%
%     subjCode = subjCodes{ss};
%
%     % run glm
%     unix(['selxavg3-sess -s ' subjCode ' -d ' data_dir ' -analysis ' analysis_name_lh_nofix ...
%         ' -no-preproc -overwrite'])
%     unix(['selxavg3-sess -s ' subjCode ' -d ' data_dir ' -analysis ' analysis_name_rh_nofix ...
%         ' -no-preproc -overwrite'])
%
% end
