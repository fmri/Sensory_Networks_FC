%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to calculate the surface area of each ROI
% in subject space using mris_anatomical_stats from freesurfer
% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;


%% Initialize Key Variables 
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
ROI_names = {'tgPCS', 'FO', 'CO', 'cIFSG', 'cmSFG', 'pAud', ...
        'sPCS', 'iPCS', 'midIFS', 'pVis', ...
        'sm_aINS', 'sm_preSMA', 'sm_dACC', 'sm_sPCS', 'sm_iPCS', 'sm_midFSG'};

N_subj = length(subjCodes);
N_ROIs = length(ROI_names);
ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
outdir = [ROI_dir 'ROI_surfarea_subjspace/'];
hemis = {'lh', 'rh'};

%% Loop over subjs, ROIs, hemis and calc surface areas
% Already ran
% parfor ss = 1:N_subj
%     subjCode = subjCodes{ss};
%     for rr = 1:N_ROIs
%         ROI_name = ROI_names{rr};
%         for hh = 1:2
%             hemi = hemis{hh};
% 
%             file_base = [ROI_dir subjCode '_' ROI_name '_' hemi];
%             outfile = [outdir subjCode '_' ROI_name '_' hemi '_indiv_anatstats.tsv'];
%             outfile_repl = [outdir subjCode '_' ROI_name '_' hemi '_replacement_indiv_anatstats.tsv'];
%             if isfile(outfile) || isfile(outfile_repl)
%                 disp(['anatomica stat file already exists for ' subjCode ROI_name hemi ' ... skipping']);
%                 continue
%             elseif isfile([file_base '_replacement.label']) % replacement labels have not yet been transformed back to subj space, so do this now
%                 unix(['mri_label2label --srclabel ' file_base '_replacement.label' ' --srcsubject fsaverage --trglabel ' ...
%                       file_base '_replacement_indiv.label ' ' --trgsubject ' subjCode ' --regmethod surface --hemi ' hemi]) 
%                 unix(['mris_anatomical_stats -l ' file_base '_replacement_indiv.label -b -f ' outfile_repl ' ' subjCode ' ' hemi])
%             elseif isfile([file_base '_indiv.label']) % already 
%                 unix(['mris_anatomical_stats -l ' file_base '_indiv.label -b -f ' outfile ' ' subjCode ' ' hemi])
%             elseif isfile([file_base '.label']) % This label has not yet been transformed back to subjs space, so do it now
%                 unix(['mri_label2label --srclabel ' file_base '.label' ' --srcsubject fsaverage --trglabel ' ...
%                       file_base '_indiv.label --trgsubject ' subjCode ' --regmethod surface --hemi ' hemi]) 
%                 unix(['mris_anatomical_stats -l ' file_base '_indiv.label -b -f ' outfile ' ' subjCode ' ' hemi])
%             elseif isfile([ROI_dir subjCode '_' hemi '_' ROI_name '.label']) % sm ROIs have a different string pattern and have not been converted to subj space, so do this now
%                 unix(['mri_label2label --srclabel ' ROI_dir subjCode '_' hemi '_' ROI_name '.label --srcsubject fsaverage --trglabel ' ...
%                       file_base '_indiv.label --trgsubject ' subjCode ' --regmethod surface --hemi ' hemi]) 
%                 unix(['mris_anatomical_stats -l ' file_base '_indiv.label -b -f ' outfile ' ' subjCode ' ' hemi])
%             else
%                 error(['No label file found for ' subjCode ' ' ROI_name ' ' hemi]);
%             end
% 
%         end
%     end
% end

%% Take mean surface area for each ROI
ROI_stats_files = {dir(outdir).name};
SA_data = nan(N_subj, N_ROIs, 2);

count = 0;

for rr = 1:N_ROIs
    ROI_name = ROI_names{rr};
    for hh = 1:2
        hemi = hemis{hh};
        if contains(ROI_name, 'sm_')
            ROI_files_curr = ROI_stats_files(contains(ROI_stats_files, ROI_name) & contains(ROI_stats_files, hemi));
        else
            ROI_files_curr = ROI_stats_files(contains(ROI_stats_files, ROI_name) & contains(ROI_stats_files, hemi) & ~contains(ROI_stats_files, '_sm_'));
        end
        if length(ROI_files_curr)~=N_subj
            keyboard
        end
        assert(length(ROI_files_curr)==N_subj);
        for ff = 1:length(ROI_files_curr)
            data = readlines([outdir ROI_files_curr{ff}]);
            data_clean = strsplit(data(61), ' ');
            if str2double(data_clean(2))<50
                count = count + 1;
            end
            SA_data(ff,rr,hh) = str2double(data_clean(3));
        end
    end
end

%% 
means_per_hemi = table(squeeze(mean(SA_data, [1])), 'RowNames', ROI_names)
means = squeeze(mean(SA_data, [1,3]));
SDs_per_hemi = table(squeeze(std(SA_data,0,1)), 'RowNames', ROI_names)
SDs = squeeze(std(SA_data, 0, [1,3]));


