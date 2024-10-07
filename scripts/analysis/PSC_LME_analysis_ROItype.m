%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to perform a linear fixed effects model for
% percent signal change for the different task contrasts
% Tom Possidente - September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/');
ccc;

%% Load subject info
load('/projectnb/somerslab/tom/projects/spacetime_network/data/behavioral/behavioral/behavioral_percent_correct_data.mat', 'perc_correct_all', 'subjCodes')
perc_correct_all = perc_correct_all(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}),1:4);
subjCodes_pcorrect = subjCodes(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}));

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.('spacetimeRuns'),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}));
N = length(subjCodes);

%% Load PSC data
load('/projectnb/somerslab/tom/projects/spacetime_network/data/PSC_results_perrun.mat', 'psc_results', 'ROIs', 'contrasts', 'subjCodes');

%% Prepare LME table
hemis = {'lh','rh'};
modalities  = {'visual', 'visual', 'auditory', 'auditory'};
domains = {'spatial', 'temporal', 'spatial', 'temporal'};
%ROIs = {'aINS', 'preSMA', 'ppreCun', 'dACC', ... % multisensory
%    'sPCS', 'iPCS', 'midIFS', 'pVis' ... % visual
%    'tgPCS', 'cIFSG', 'pAud', 'CO', 'FO', 'cmSFG'}; % auditory
ROIs_visual = {'sPCS', 'iPCS', 'midIFS'};
ROIs_auditory = {'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG'};
runs = {'run1', 'run2', 'run3', 'run4', 'run5', 'run6'};
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
            prc_correct = perc_correct_all{ss,cc};
            for rr = 1:dims(4)
                ROI = ROIs{rr};
                if ismember(ROI, ROIs_visual)
                    ROItype = 'visualROI';
                elseif ismember(ROI, ROIs_auditory)
                    ROItype = 'auditoryROI';
                else
                    continue;
                end
                for run = 1:dims(5)
                    run_text = runs{run};
                    psc = psc_results(ss,hh,cc,rr,run);
                    if isnan(psc)
                        continue
                    else
                        data_table = [data_table; {psc, ss, hemi, ROItype, modality, domain, run}];
                    end
                end
            end
        end
    end
end

data_table.Properties.VariableNames = {'PSC', 'subject', 'hemisphere', 'ROItype', 'modality', 'domain', 'run'};
data_table.hemisphere = categorical(data_table.hemisphere);
data_table.ROItype = categorical(data_table.ROItype);
data_table.modality = categorical(data_table.modality);
data_table.domain = categorical(data_table.domain);
data_table.run = categorical(data_table.run);


%% Fit LME

lme = fitlme(data_table, ['PSC ~ 1 + modality * domain * ROItype + (1 + ROItype + modality + domain | subject) ' ...
    ' + (1 + ROItype + modality + domain | hemisphere)'])


% Do we expect different subjects to have different PSC intercepts/slopes
% based on ROItype, modality,  and domain?
% - Yes to all?

% Do we expect different hemispheres to have different PSC intercepts/slopes
% based on ROItype, modality,  and domain?
% - yes modality, maybe domain, maybe ROItype

% Do we expect different runs to have different PSC intercepts/slopes
% based on ROItype, modality,  and domain?
% - No, just use single intercept


%% Permutation testing - generate null distributions
original_tstats = lme.Coefficients.tStat;
tstat_ROItype_modality = original_tstats(5);
tstat_ROItype_domain = original_tstats(6);
tstat_ROItype_modality_domain = original_tstats(8);

original_coef_names = lme.Coefficients.Name;
iterations = 5000;

null_dist_tstats = nan(iterations, length(original_tstats([5,6,8])));

tic;
parfor ii = 1:iterations
    permuted_table = table();

    % Shuffle modality, domain, and ROI type labels for each subject
    for ss = 1:dims(1)
        subj_data = data_table(data_table.subject == ss,:);
        num_obs = height(subj_data);
        perm_inds = [randperm(num_obs); randperm(num_obs); randperm(num_obs)]';
        subj_data.ROItype = subj_data.ROItype(perm_inds(:,1));
        subj_data.modality = subj_data.modality(perm_inds(:,2));
        subj_data.domain = subj_data.domain(perm_inds(:,3));
        permuted_table = [permuted_table; subj_data];
    end

    perm_lme = fitlme(permuted_table, ['PSC ~ 1 + modality * domain * ROItype + (1 + ROItype + modality + domain | subject) ' ...
    ' + (1 + ROItype + modality + domain | hemisphere)']);
    null_dist_tstats(ii,:) = perm_lme.Coefficients.tStat([5,6,8]);
end

null_dist_tstats_max_pos = max(null_dist_tstats,[],2);
null_dist_tstats_max_neg = min(null_dist_tstats,[],2);
null_dist_tstats_max_abs = max(abs(null_dist_tstats),[],2);

pval_ROItype_modality = (sum(tstat_ROItype_modality < null_dist_tstats_max_pos)+1) / iterations
pval_ROItype_domain = (sum(tstat_ROItype_domain > null_dist_tstats_max_neg)+1) / iterations
pval_ROItype_modality_domain_ = (sum(tstat_ROItype_modality_domain < null_dist_tstats_max_neg)+1) / iterations

toc


