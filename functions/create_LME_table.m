function LME_table = create_LME_table(psc_data, prc_correct_data, hemis, modalities, domains, ROItype_names, make_categorical, bad_subj_mask_percond)
%CREATE_LME_TABLE 
%
% Inputs:
%   psc_data: 4D double sxhxcxr where s is subjects, h is hemispheres, c is
%             conditions, and r is ROI types
%   prc_correct_data: 2D double sxc where s is subjects and c is conditions
%   hemis: cell array - ordered hemisphere names (ex. {'lh','rh'})
%   modalities: cell array - ordered modality names for each condition
%   domains: cell array - ordered domains for each condition
%   ROItype_names: cell array - ordered ROI type names for each ROItype 
%   make_categorical: logical - make the hemisphere, modality, domain, and ROItype
%   bad_subj_inds_pertask: 2D logical sxc where s is subjects and c is conditions
%
% Outputs:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Inputs
assert(nargin >= 5, 'Need 4 or more input arguments to creat_LME_table funcion');

%% Loop through data and create LME table
dims = size(psc_data);
LME_table = table();
for ss = 1:dims(1)
    for hh = 1:dims(2)
        hemi = hemis{hh};
        for cc = 1:dims(3)
            if bad_subj_mask_percond(ss,cc)
                continue;
            end
            modality = modalities{cc};
            domain = domains{cc};
            prc_correct = prc_correct_data(ss,cc);
            for rr = 1:dims(4)
                ROI = ROItype_names{rr};
                psc = psc_data(ss,hh,cc,rr);
                if isnan(psc)
                    continue
                else
                    LME_table = [LME_table; {psc, ss, hemi, ROI, modality, domain, prc_correct}];
                end
            end
        end
    end
end

LME_table.Properties.VariableNames = {'PSC', 'subject', 'hemisphere', 'ROItype', 'modality', 'domain', 'perc_correct'};

if make_categorical
    LME_table.hemisphere = categorical(LME_table.hemisphere);
    LME_table.ROItype = categorical(LME_table.ROItype);
    LME_table.modality = categorical(LME_table.modality);
    LME_table.domain = categorical(LME_table.domain);
end

end

