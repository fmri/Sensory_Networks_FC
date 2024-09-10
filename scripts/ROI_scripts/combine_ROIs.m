%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to load label files and combine them into
%%% one larger label file (useful for taking all posterior visual ROIs and
%%% combining them into a single pVis)
%%% Tom Possidente - August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Get Subj Codes

projectDir = '/projectnb/somerslab/tom/projects/spacetime_network/';
experiment_name = 'spacetime';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'SL', 'AH', 'RR'})); % rejected subjs
N = length(subjCodes);

%% Set key variables
pVis_ROIs = {'VOT', 'LOT', 'aIPS', 'pIPS'};
%pVis_ROIs = {'VOT', 'LOT', 'DO', 'aIPS', 'pIPS'};

ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
ROI_files = {dir(ROI_dir).name};
hemis = {'lh', 'rh'};

%% Loop through subjs and combine posterior visual ROIs 

for ss = 1:N
    subjCode = subjCodes{ss};
    for hh = 1:2
        hemi = hemis{hh};
        subj_pvis_ROI_fnames = ROI_files(contains(ROI_files, pVis_ROIs) & ~contains(ROI_files, '_indiv') ...
                                         & contains(ROI_files, ['_' hemi]) & contains(ROI_files, [subjCode '_']));
        if strcmp(subjCode, 'NS') % due to "NS" being in the label aINS, we have to take out the wrong aINS for this subj
            bad_aINS = contains(subj_pvis_ROI_fnames, 'aINS');
            subj_pvis_ROI_fnames = subj_pvis_ROI_fnames(~bad_aINS);
        end

        replacement_ROIs = contains(subj_pvis_ROI_fnames, 'replacement.label');
        if any(replacement_ROIs)
            assert(sum(replacement_ROIs)<=3, ['more than one replacement ROI found in pVis ROIs for subj ' subjCode])
            repl_ROI_name = strsplit(subj_pvis_ROI_fnames{replacement_ROIs}, '_');
            repl_ROI_name = repl_ROI_name{2};
            bad_nonrepl = contains(subj_pvis_ROI_fnames, '.nii') | contains(subj_pvis_ROI_fnames, [repl_ROI_name '_' hemi '.label']);
            subj_pvis_ROI_fnames = subj_pvis_ROI_fnames(~bad_nonrepl);
        end

        assert(length(subj_pvis_ROI_fnames)==length(pVis_ROIs), ['Wrong number of pVis ROIs found for subj ' subjCode ...
                                                                 '. Expected ' num2str(length(pVis_ROIs)) ' found ' ...
                                                                 num2str(length(subj_pvis_ROI_fnames))]);
        % Loop through ROIs and combine vertices
        pVis = [];
        for rr = 1:length(pVis_ROIs)
            ROI = readtable([ROI_dir subj_pvis_ROI_fnames{rr}], "FileType","text");
            pVis = [pVis; table2array(ROI)];
        end
        
        % Save out pVis ROI
        label_fname = [ROI_dir subjCode '_pVis_mod_' hemi '.label'];
        [ROI_unique,~,~] = unique(pVis,'rows');
        label_file = fopen(label_fname,'w');
        fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(ROI_unique,1)) '\n']);
        writematrix(ROI_unique, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
        fclose(label_file);

    end
end














