%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to perform a linear fixed effects model for
% percent signal change for the different task contrasts
% Tom Possidente - September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load subject info
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.('spacetimeRuns'),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}));
N = length(subjCodes);

%% Load PSC data
load('/projectnb/somerslab/tom/projects/spacetime_network/data/PSC_results.mat', 'psc_results', 'ROIs', 'contrasts', 'subjCodes');

%% Prepare LME table
hemis = {'lh','rh'};
modalities  = {'visual', 'visual', 'auditory', 'auditory'};
domains = {'spatial', 'temporal', 'spatial', 'temporal'};
%ROIs = {'aINS', 'preSMA', 'ppreCun', 'dACC', ... % multisensory
%    'sPCS', 'iPCS', 'midIFS', 'pVis' ... % visual
%    'tgPCS', 'cIFSG', 'pAud', 'CO', 'FO', 'cmSFG'}; % auditory
ROIs_visual = {'sPCS', 'iPCS', 'midIFS'};
ROIs_auditory = {'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG'};
ROIs_use = {'sPCS', 'iPCS', 'midIFS' ... % visual
       'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG'}; % auditory
dims = size(psc_results);

data_table = table();
for ss = 1:dims(1)
    for hh = 1:dims(2)
        hemi = hemis{hh};
        for cc = 1:dims(3)
            if cc > 4 % only using 4 tasks: sV, tV, sA, tA
                continue
            end
            modality = modalities{cc};
            domain = domains{cc};
            for rr = 1:dims(4)
                ROI = ROIs{rr};
                if ismember(ROI, ROIs_visual)
                    ROItype = 'visualROI';
                elseif ismember(ROI, ROIs_auditory)
                    ROItype = 'auditoryROI';
                else
                    continue;
                end

                psc = psc_results(ss,hh,cc,rr);
                if isnan(psc)
                    continue
                else
                    data_table = [data_table; {psc, ss, hemi, ROItype, modality, domain}];
                end
            end
        end
    end
end

data_table.Properties.VariableNames = {'PSC', 'subject', 'hemisphere', 'ROItype', 'modality', 'domain'};
data_table.hemisphere = nominal(data_table.hemisphere);
data_table.ROItype = nominal(data_table.ROItype);
data_table.modality = nominal(data_table.modality);
data_table.domain = nominal(data_table.domain);

%% Fit LME

lme = fitlme(data_table, ['PSC ~ 1 + modality * domain * ROItype + (1 + modality + domain + ROItype | subject) ' ...
    ' + (1 + modality + domain + ROItype | hemisphere)'])









