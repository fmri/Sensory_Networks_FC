%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to take the full HCP Glasser atlas ROIs
%%% annotation files and create a new annotation file containing just the
%%% relevant domain general parcels.
%%% Tom Possidente - April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[verts_ref, labels_ref, ctable_ref] = read_annotation('lh.HCP-MMP1.annot'); % Read in lh.aparc.annot file for reference annotation file structure
keep = {'L_s6-8_ROI', 'L_i6-8_ROI', 'L_8C_ROI', 'L_IFJp_ROI', 'L_p9-46v_ROI', 'L_6r_ROI', 'L_p10p_ROI', 'L_a10p_ROI', 'L_11l_ROI',...
        'L_a47r_ROI', 'L_p47r_ROI', 'L_a9-46v_ROI', 'L_FOP5_ROI', 'L_AVI_ROI', 'L_TE1m_ROI', 'L_TE1p_ROI', 'L_AIP_ROI', 'L_IP2_ROI',...
        'L_IP1_ROI', 'L_LIPd_ROI', 'L_MIP_ROI', 'L_PGs_ROI', 'L_PFm_ROI', 'L_POS2_ROI', 'L_SCEF_ROI', 'L_8BM_ROI', 'L_a32pr_ROI', 'L_d32_ROI'};
keep_mask = ismember(ctable_ref.struct_names,keep);

ctable_ref.numEntries = length(keep);
ctable_ref.struct_names = ctable_ref.struct_names(keep_mask);
ctable_ref.table = ctable_ref.table(keep_mask,:);

write_annotation('lh.HCP-MMP1_assemMD.annot', verts_ref, labels_ref, ctable_ref);

%%
[verts_ref, labels_ref, ctable_ref] = read_annotation('rh.HCP-MMP1.annot'); % Read in lh.aparc.annot file for reference annotation file structure
keep = {'R_s6-8_ROI', 'R_i6-8_ROI', 'R_8C_ROI', 'R_IFJp_ROI', 'R_p9-46v_ROI', 'R_6r_ROI', 'R_p10p_ROI', 'R_a10p_ROI', 'R_11l_ROI',...
        'R_a47r_ROI', 'R_p47r_ROI', 'R_a9-46v_ROI', 'R_FOP5_ROI', 'R_AVI_ROI', 'R_TE1m_ROI', 'R_TE1p_ROI', 'R_AIP_ROI', 'R_IP2_ROI',...
        'R_IP1_ROI', 'R_LIPd_ROI', 'R_MIP_ROI', 'R_PGs_ROI', 'R_PFm_ROI', 'R_POS2_ROI', 'R_SCEF_ROI', 'R_8BM_ROI', 'R_a32pr_ROI', 'R_d32_ROI'};
keep_mask = ismember(ctable_ref.struct_names,keep);

ctable_ref.numEntries = length(keep);
ctable_ref.struct_names = ctable_ref.struct_names(keep_mask);
ctable_ref.table = ctable_ref.table(keep_mask,:);

write_annotation('rh.HCP-MMP1_assemMD.annot', verts_ref, labels_ref, ctable_ref);