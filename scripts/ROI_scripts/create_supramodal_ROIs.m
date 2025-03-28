%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create individual subject supramodal
%%% ROIs using the intersection of WM aud and vis vertices within
%%% particular search spaces defined from probabilistic maps. Using
%%% localizer data.
%%%
%%% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load Subj Codes
reject_subjs = {'RR', 'AH', 'SL'};
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.('spacetimeRuns'),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, reject_subjs));
n_subj = length(subjCodes);
badlist = {};

%% Set key variables
ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
subj_data_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/unpacked_data_nii_fs_localizer/';
hemis = {'lh', 'rh'};

t_thresh = 2;
t_thresh_2 = 2;

supramodal_ROIs = {'sPCS', 'iPCS', 'midFSG', 'aINS', 'preSMA', 'dACC'} %'aIPS', 'pIPS'};
n_ROIs = length(supramodal_ROIs);

probROI_data = cell(n_ROIs, 2);

%% Load probabilistic ROIs (search spaces for individual ROIs)
for rr = 1:n_ROIs
    for hh = 1:2
        path = [ROI_dir hemis{hh} '_intersect_' supramodal_ROIs{rr} '_supra_vA-vP+aA-aP_probabilistic_thresh5_dilate5.label'];
        probROI_data{rr,hh} = readtable(path, 'FileType','text');
    end
end

%% Loop through subjs and create ROIs
num_imperfect_ROIs = 0;
num_0_ROIs = 0;
num_bad_ROIs = 0;
lengths = nan(n_subj, n_ROIs, 2);
for ss = 1:n_subj
    for hh = 1:2
        % Load subj intersection map
        path_vA_vP = [subj_data_dir subjCodes{ss} '/localizer/localizer_contrasts_' hemis{hh} '/vA-vP/t.nii.gz'];
        path_aA_aP = [subj_data_dir subjCodes{ss} '/localizer/localizer_contrasts_' hemis{hh} '/aA-aP/t.nii.gz'];
        path_V_A = [subj_data_dir subjCodes{ss} '/localizer/localizer_contrasts_' hemis{hh} '/V-A/t.nii.gz'];
        vA_vP = MRIread(path_vA_vP);
        aA_aP = MRIread(path_aA_aP);
        V_A = MRIread(path_V_A);
        for rr = 1:n_ROIs
            probROI = probROI_data{rr,hh};
            vA_vP_inROI = vA_vP.vol(probROI.Var1+1); % add one to get label inds to line up with nii inds
            aA_aP_inROI = aA_aP.vol(probROI.Var1+1); % add one to get label inds to line up with nii inds
            V_A_inROI = V_A.vol(probROI.Var1+1); % add one to get label inds to line up with nii inds

            intersection_verts = (vA_vP_inROI > t_thresh) & (vA_vP_inROI > t_thresh) & (abs(V_A_inROI) < t_thresh_2);

            if sum(sum(intersection_verts)<100)
                if sum(intersection_verts)==0
                    num_0_ROIs = num_0_ROIs + 1;
                end
                num_imperfect_ROIs = num_imperfect_ROIs + 1;
                badlist{end+1} = [subjCodes{ss} '_sm_' supramodal_ROIs{rr} '_' hemis{hh}];
                thresh_alt = t_thresh;
                thresh_alt_2 = t_thresh_2;
                while sum(intersection_verts)<100
                    if sum((vA_vP_inROI > thresh_alt) & (vA_vP_inROI > thresh_alt)) < sum(abs(V_A_inROI) < thresh_alt_2)
                        thresh_alt = thresh_alt - 0.01;
                    else
                        thresh_alt_2 = thresh_alt_2 + 0.01;
                    end
                    intersection_verts = (vA_vP_inROI > thresh_alt) & (vA_vP_inROI > thresh_alt) & (abs(V_A_inROI) < thresh_alt_2);
                end
                disp([subjCodes{ss} ' ' hemis{hh} ' ' supramodal_ROIs{rr} ' needed a t=' num2str(thresh_alt) ', t2= ' num2str(thresh_alt_2) ' threshold to produce a 100 vertex ROI'])
                if any(vA_vP_inROI(intersection_verts)<0) | any(vA_vP_inROI(intersection_verts)<0)
                    num_bad_ROIs = num_bad_ROIs + 1;
                    disp([num2str(sum(vA_vP_inROI(intersection_verts)<0)) ' vA-vP vertices have a t-stat below 0']);
                    disp([num2str(sum(aA_aP_inROI(intersection_verts)<0)) ' aA-aP vertices have a t-stat below 0']);
                end
            end

            intersection_verts_inds = find(intersection_verts); % minus one to convert back to label indices
            label_table = probROI(intersection_verts_inds,:);
            
            % Make label file
            % label_fname = [ROI_dir subjCodes{ss} '_' hemis{hh} '_sm_' supramodal_ROIs{rr} '2.label'];
            % label_file = fopen(label_fname,'w');
            % fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label_table,1)) '\n']);
            % writematrix(table2array(label_table), label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
            % fclose(label_file);

        end

    end
end






