function LME_table_out = vis_aud_2contrast(LME_table)
%VIS_AUD_2CONTRAST  create double contrast of visWM-visPassive -
%audWM-audPassive from LME table of PSCs. 

LME_table_out = table();

subjs = unique(LME_table.subject);
hemis = {'lh','rh'};
ROIs = unique(LME_table.ROItype);
domains = unique(LME_table.domain);

for ss = 1:length(subjs)
    subj_mask = ismember(LME_table.subject, subjs(ss));
    for rr = 1:length(ROIs)
        ROI_mask = ismember(LME_table.ROItype, ROIs{rr});
        for dd = 1:length(domains)
            domain_mask = ismember(LME_table.domain, domains{dd});
            for hh = 1:length(hemis)
                hemis_mask = ismember(LME_table.hemisphere, hemis{hh});
                full_mask = subj_mask & ROI_mask & domain_mask & hemis_mask;
                va_rows = LME_table(full_mask,:);
                if height(va_rows)~=2
                    disp(['rows in vis-aud contrast table not equal to 2 for subj ' num2str(ss), ' ROI ' ROIs{rr}, ' ' hemis{hh}]);
                    continue;
                end
                new_PSC = va_rows.PSC(ismember(va_rows.modality,'visual')) - va_rows.PSC(ismember(va_rows.modality,'auditory'));
                assert(~isnan(new_PSC), 'new PSC is nan');
                new_pcorrect = va_rows.perc_correct(ismember(va_rows.modality,'visual')) - va_rows.perc_correct(ismember(va_rows.modality,'auditory'));
                if abs(new_pcorrect) > 1
                    keyboard;
                end
                new_row = va_rows(1,:);
                new_row.PSC = new_PSC;
                new_row.perc_correct = new_pcorrect;
                new_row.modality = 'v-a';
                LME_table_out = [LME_table_out; new_row];
            end
        end
    end
end

LME_table_out.modality = cellstr(LME_table_out.modality); % convert so it's in the right format to make this a categorical array later

if height(LME_table_out) ~= height(LME_table)/2
    disp('Rows in LME_table_out are not half of rows in LME_table. Likely caused by missing visual or auditory conditions for some subjs/ROIs');
end



end

