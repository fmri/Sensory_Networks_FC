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
%pVis_ROIs = {'VOT', 'LOT', 'aIPS', 'pIPS'};
%ROIs_combine = {'sPCS', 'iPCS', 'midIFS'};
%ROIs_combine = {'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG'};
ROIs_combine = {'aINS', 'preSMA', 'dACC'};
%combined_name = 'vbiased_frontal';
%combined_name = 'abiased_frontal';
combined_name = 'MD_frontal';

ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
ROI_files = {dir(ROI_dir).name};
hemis = {'lh', 'rh'};

load('/projectnb/somerslab/tom/projects/spacetime_network/data/missing_ROIs.mat', 'missing_ROIs')

%% Loop through subjs and combine posterior visual ROIs 

for ss = 1:N
    subjCode = subjCodes{ss};
    for hh = 1:2
        hemi = hemis{hh};
        combine_ROIs_fnames = ROI_files(contains(ROI_files, ROIs_combine) & ~contains(ROI_files, '_indiv') ...
                                         & contains(ROI_files, ['_' hemi]) & contains(ROI_files, [subjCode '_']));
        if strcmp(subjCode, 'NS') % due to "NS" being in the label aINS, we have to take out the wrong aINS for this subj
            bad_aINS = contains(combine_ROIs_fnames, 'aINS') & ~contains(combine_ROIs_fnames, 'NS_aINS');
            combine_ROIs_fnames = combine_ROIs_fnames(~bad_aINS);
        end

        replacement_ROI_inds = contains(combine_ROIs_fnames, 'replacement.label');
        if any(replacement_ROI_inds)
            replacement_ROIs = combine_ROIs_fnames(replacement_ROI_inds);
            for rr = 1:length(replacement_ROIs)
                repl_ROI_name = strsplit(replacement_ROIs{rr}, '_');
                repl_ROI_name = repl_ROI_name{2};
                bad_nonrepl = contains(combine_ROIs_fnames, '.nii') | contains(combine_ROIs_fnames, [repl_ROI_name '_' hemi '.label']);
                combine_ROIs_fnames = combine_ROIs_fnames(~bad_nonrepl);
            end
        end
        
        if length(combine_ROIs_fnames)~=length(ROIs_combine) % make sure correct number of ROIs found
            num_missing = sum( contains(missing_ROIs, subjCode) & contains(missing_ROIs, hemi) & contains(missing_ROIs, ROIs_combine) );
            assert(length(combine_ROIs_fnames)==length(ROIs_combine)-num_missing, ['Wrong number of ROIs found for subj ' subjCode ...
                                                                 '. Expected ' num2str(length(ROIs_combine)) ' found ' ...
                                                                 num2str(length(combine_ROIs_fnames))]);
        else
            num_missing = 0;
        end

        % Loop through ROIs and combine vertices
        combined = [];
        for rr = 1:length(ROIs_combine)-num_missing
            ROI = readtable([ROI_dir combine_ROIs_fnames{rr}], "FileType","text");
            combined = [combined; table2array(ROI)];
        end
        
        % Save out combined ROI
        label_fname = [ROI_dir subjCode '_' combined_name '_' hemi '.label'];
        [ROI_unique,~,~] = unique(combined,'rows');
        label_file = fopen(label_fname,'w');
        fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(ROI_unique,1)) '\n']);
        writematrix(ROI_unique, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
        fclose(label_file);

    end
end














