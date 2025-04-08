%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to calculate the surface area of each ROI
% in subject space using mris_anatomical_stats from freesurfer
% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/'));
ccc;


%% Initialize Key Variables
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
ROI_names = {'tgPCS', 'FO', 'CO', 'cIFSG', 'cmSFG', 'pAud', ...
    'sPCS', 'iPCS', 'midIFS', 'pVis', ...
    'sm_aINS', 'sm_preSMA', 'sm_dACC', 'sm_sPCS', 'sm_iPCS', 'sm_midFSG'};

N_subj = length(subjCodes);
N_ROIs = length(ROI_names);
ROI_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs/';
ROI_final_fs_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs_final/labels_fsaverage/';
ROI_final_subj_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs_final/labels_subjspace/';
outdir = [ROI_dir 'ROI_surfarea_subjspace/'];
hemis = {'lh', 'rh'};

%% Loop over subjs, ROIs, hemis and calc surface areas
for ss = 1:N_subj
    subjCode = subjCodes{ss};
    for hh = 1:2
        hemi = hemis{hh};

        file_path = [ROI_dir hemi '.' subjCode '_avsm_ROIs.annot'];
        
        if ~isfolder([ROI_final_fs_dir subjCode])
            unix(['mkdir ' ROI_final_fs_dir subjCode])
        end
        if ~isfolder([ROI_final_subj_dir subjCode])
            unix(['mkdir ' ROI_final_subj_dir subjCode])
        end


        unix(['mri_annotation2label --subject fsaverage --hemi ' hemi ' --outdir ' ROI_final_fs_dir subjCode '/ --annotation ' file_path]);

        for rr = 1:N_ROIs
            ROI_name = ROI_names{rr};
            outfile = [outdir subjCode '_' ROI_name '_' hemi '_anatstats.tsv'];
            unix(['mri_label2label --srclabel ' ROI_final_fs_dir subjCode '/' hemi '.' ROI_name '.label ' '--srcsubject fsaverage --trglabel ' ...
                ROI_final_subj_dir subjCode '/' hemi '.' ROI_name '.label ' ' --trgsubject ' subjCode ' --regmethod surface --hemi ' hemi]) 
            unix(['mris_anatomical_stats -l ' ROI_final_subj_dir subjCode '/' hemi '.' ROI_name '.label ' ' -b -f ' outfile ' ' subjCode ' ' hemi])
        end

    end
end


%% Take mean surface area for each ROI
ROI_stats_files = {dir(outdir).name};
SA_data = nan(N_subj, N_ROIs, 2);

count = 0;

for rr = 1:N_ROIs
    ROI_name = ROI_names{rr};
    for hh = 1:2
        hemi = hemis{hh};
        if contains(ROI_name, 'sm_')
            ROI_files_curr = ROI_stats_files(contains(ROI_stats_files, ROI_name) & contains(ROI_stats_files, hemi) & ~contains(ROI_stats_files, 'indiv'));
        else
            ROI_files_curr = ROI_stats_files(contains(ROI_stats_files, ROI_name) & contains(ROI_stats_files, hemi) & ~contains(ROI_stats_files, '_sm_') &  ~contains(ROI_stats_files, 'indiv'));
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

%% Get means and SDs
means_per_hemi = table(squeeze(mean(SA_data, [1])), 'RowNames', ROI_names)
means = squeeze(mean(SA_data, [1,3]));
SDs_per_hemi = table(squeeze(std(SA_data,0,1)), 'RowNames', ROI_names)
SDs = squeeze(std(SA_data, 0, [1,3]));


