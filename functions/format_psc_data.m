function [psc_results, perc_correct_all, hemis, modality_out, domain_out, ROIs, bad_subjs] = ...
                    format_psc_data(modalities_use, domains_use, use_ROItypes, include_MD, localizer_data)
%FORMAT_PSC_DATA 
%
%

%% Load percent correct data
if ~localizer_data
    load('/projectnb/somerslab/tom/projects/spacetime_network/data/behavioral/behavioral/behavioral_percent_correct_data.mat', 'perc_correct_all', 'subjCodes')
    perc_correct_all = table2array(perc_correct_all(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}),1:4)); % cut out bad subjs and tactile % correct data
    perc_correct_all(:,5:6) = nan(size(perc_correct_all,1),2); % add 2 columns of nans for passive conditions for which there is no task (and therefore no % correct data)
    subjCodes_pcorrect = subjCodes(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'})); % cut rejected subjs from subjCodes as well
    task_str = 'spacetime';
    psc_fname = 'PSC_results.mat';
    modality_order  = {'visual', 'visual', 'auditory', 'auditory', 'visual', 'auditory'};
    domain_order = {'spatial', 'temporal', 'spatial', 'temporal', 'passive', 'passive'};
else
    task_str = 'x1WayLocalizer';
    psc_fname = 'PSC_results_localizer.mat';
    modality_order  = {'visual','auditory', 'auditory', 'visual'};
    domain_order = {'active', 'active', 'passive', 'passive'};
end

%% Load subject info
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([task_str 'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'AH', 'SL', 'RR', 'AI'}));
N = length(subjCodes);
if ~localizer_data 
    assert(all(isequal(subjCodes,subjCodes_pcorrect))); % make sure same subjCodes in same order
end

%% Load PSC data
subjCodes_pre = subjCodes;
load(['/projectnb/somerslab/tom/projects/spacetime_network/data/PSCs/' psc_fname], 'psc_results', 'ROIs', 'contrasts', 'subjCodes');
use_ROIs = ~ismember(ROIs,{'ppreCun'});
psc_results = psc_results(:,:,:,use_ROIs);
ROIs = ROIs(use_ROIs);
assert(all(isequal(subjCodes,subjCodes_pre))); % make sure same subjCodes in same order

N_conds = size(psc_results, 3);
N_ROIs = size(psc_results, 4);

%% Get conditions being used
modality_mask = ismember(modality_order, modalities_use);
domain_mask = ismember(domain_order, domains_use);
c_mask = modality_mask & domain_mask;

%% Create badsubjs array
bad_subjs = false(N,N_conds);
if ~localizer_data
    badsubjs_auditory_tasks = {'LN', 'GG', 'TP'};
    badsubjs_visual_tasks = {'LA', 'RT'};
    for cc = 1:N_conds
        if strcmp('passive', domain_order{cc})
            continue;
        elseif strcmp('visual', modality_order{cc})
            bad_subjs(:,cc) = ismember(subjCodes, badsubjs_visual_tasks);
        elseif strcmp('auditory', modality_order{cc})
            bad_subjs(:,cc) = ismember(subjCodes, badsubjs_auditory_tasks);
        else
            error('unrecognized condition');
        end
    end
end

%% Get ROIs being used
ROIs_MD = {'aINS', 'dACC', 'preSMA'};
if include_MD
    ROI_mask = true(N_ROIs,1);
else
    ROI_mask = ~ismember(ROIs, ROIs_MD);
    ROIs = ROIs(ROI_mask);
end

%% Cut psc_results down to conditions and ROIs being used
psc_results = psc_results(:,:,c_mask,ROI_mask);

%% Configure ROItypes being used
if use_ROItypes
    ROIs_visual = {'sPCS', 'iPCS', 'midIFS'};
    ROIs_auditory = {'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG'};
    ROIs_MD = {'aINS', 'dACC', 'preSMA'};
    
    ROIs(ismember(ROIs, ROIs_visual)) = repmat({'visualROI'}, length(ROIs_visual), 1);
    ROIs(ismember(ROIs, ROIs_auditory)) = repmat({'auditoryROI'}, length(ROIs_auditory), 1);
    if include_MD
        ROIs(ismember(ROIs, ROIs_MD)) = repmat({'MDROI'}, length(ROIs_MD), 1);
    end
else
    ROItype_names = ROIs;
end

%% Create LME table
hemis = {'lh','rh'};
if ~localizer_data
    perc_correct_all = perc_correct_all(:,c_mask); % extract % correct data from only conditions being used
else
    perc_correct_all = nan(N,N_conds-sum(~c_mask));
end
bad_subjs = bad_subjs(:,c_mask); % extract bad subj indicators only from conditions being used

modality_out = modality_order(c_mask);
domain_out = domain_order(c_mask);

end

