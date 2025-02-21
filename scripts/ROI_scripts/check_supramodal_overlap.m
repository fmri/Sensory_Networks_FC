%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to assess how much overlap there is between
% the individual subject supramodal ROIs and the individual abias/vbias ROIs
%
% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Initialize parameters
ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'}; % first 20 are localizer analysis subjs, AI is in rs only
N = length(subjCodes);
av_biased_ROIs = {'sPCS', 'iPCS', 'midIFS', 'tgPCS', 'cIFSG', 'CO', 'FO'};
supramodal_ROIs = {'sm_sPCS', 'sm_iPCS', 'sm_midFSG'};
hemis = {'lh', 'rh'};
overlaps = zeros([N, length(av_biased_ROIs), length(supramodal_ROIs), 2]);
lower_100_count = 0;
total_count = 0;
overlap_count = 0;

%% Loop through subjs and ROIs to assess overlap
for ss = 1:N
    subjCode = subjCodes{ss};
    for hh = 1:2
        hemi = hemis{hh};

        % Load annotation files for both ROI sets for this subj
        [verts_ref, labels_orig, ctable_orig] = read_annotation([ROI_dir hemi '.' subjCode '_ROIs.annot']); 
        [~        , labels_sm, ctable_sm] = read_annotation([ROI_dir hemi '.' subjCode '_smROIs2.annot']); 

        for rb = 1:length(av_biased_ROIs)
            av_ROI = av_biased_ROIs{rb};
            av_ROI_mask = labels_orig == ctable_orig.table(strcmp(ctable_orig.struct_names,av_ROI),5);
            for rs = 1:length(supramodal_ROIs)
                if rs < rb
                    continue
                end
                total_count = total_count + 1;
                sm_ROI = supramodal_ROIs{rs};
                sm_ROI_mask = labels_sm == ctable_orig.table(strcmp(ctable_sm.struct_names,sm_ROI),5);
                overlap = sum(sm_ROI_mask & av_ROI_mask);
                if overlap ~= 0
                    overlap_count = overlap_count + 1;
                    disp([subjCode ' ' hemi ' ' av_ROI ' and ' sm_ROI ' overlap ' num2str(overlap) ' vertices']);
                    if sum(sm_ROI_mask)-overlap < 100 && sum(av_ROI_mask)-overlap < 100
                        lower_100_count = lower_100_count + 1;
                    end
                end
                overlaps(ss,rb,rs,hh) = overlap;
            end
        end
    end
end






