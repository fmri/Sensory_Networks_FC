%%
function freeview_screenshots(subjIDs, contrast_list, lh_label_list, rh_label_list, colortable, use_fsaverage, ...
    save_dir, func_path, func_folder, recon_dir, analysis_name, stat_file, overlayThreshold, label_opacity, ss_suffix, ss_on, overlay_path_override)
%FREEVIEW_SCREENSHOTS - Creates command files for freeview to plot and screenshot all given labels per contrast and subj
%
% If matlab is running into errors opening freesurfer, try running this
% with freeview version 6.0 loaded.
%
% David Beeler, Tom Possidente - August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize inputs
assert(nargin >= 4, 'Not enough input arguments (subjIDs, contrast_list, lh_label_list, rh_label_list all required)');

if nargin < 5 || isempty(colortable) && (~isempty(lh_label_list) && ~isempty(rh_label_list))
    ROI = ["aINS", "preSMA", "ppreCun", "dACC", ... % multisensory
        "sPCS", "iPCS", "midIFS", "aIPS", "pIPS", "DO", "LOT", "VOT",... % visual
        "tgPCS", "cIFSG", "pAud", "CO", "FO", "cmSFG"]'; % auditory
    color = [[0 255 0]; [0 255 0]; [0 255 0]; [0 255 0]; ... % green
        [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; [0 0 255]; ... % blue
        [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]; [255 165 0]]; % orange
    colortable = table(ROI, color(:,1), color(:,2), color(:,3), 'VariableNames', {'ROI', 'c1', 'c2', 'c3'});
end

if nargin < 6 || isempty(use_fsaverage)
    use_fsaverage = true;
end

if nargin < 7 || isempty(save_dir)
    save_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/figures_images/roi_QC_screenshots/';
end

if nargin < 8 || isempty(func_path)
    func_path = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
end

if nargin < 9 || isempty(func_folder)
    func_folder = 'localizer';
end

if nargin < 10 || isempty(recon_dir)
    recon_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/';
end

if nargin < 11 || isempty(analysis_name)
    analysis_name = {'localizer_contrasts', ''};
end

if nargin < 12 || isempty(stat_file)
    stat_file = {'sig.nii.gz'};
end

if nargin < 13 || isempty(overlayThreshold)
    overlayThreshold = '1.3,3,5'; % min, mid, max
end

if nargin < 14 || isempty(label_opacity)
    label_opacity = '0.5'; % min, mid, max
end

if nargin < 15 || isempty(ss_suffix)
    ss_suffix = '';
end

if nargin < 16 || isempty(ss_on)
    ss_on = true;
end

if nargin < 17 || isempty(overlay_path_override)
    overlay_path_override = false;
end


% Set key variables
inflationType = 'inflated'; % inflated / pial
hemis = {'lh','rh'};
viewAngles = {'lateral', 'medial'}; % could add dorsal / ventral
%viewAngles = {'lateral', 'medial', 'posterior'}; % could add dorsal / ventral

labelOutline = 'false'; % true / false
annotOutline = 'true';

if overlay_path_override
    lh_analysisName = analysis_name;
    rh_analysisName = analysis_name;
else
    lh_analysisName = [analysis_name{1} '_lh_' analysis_name{2}];
    rh_analysisName = [analysis_name{1} '_rh_' analysis_name{2}];
end

if isfile([save_dir '/run_all_screenshots.sh'])
    disp('run_all_screenshots.sh already exists in the save_dir, pausing execution. Either exit the script or press continue to delete and replace the file.')
    keyboard;
    unix(['rm ' save_dir '/run_all_screenshots.sh']);
end

no_contrasts = false;
if isempty(contrast_list)
    no_contrasts = true;
    contrast_list = cell(length(subjIDs),1);
    contrast_list(:) = {''};
end

for ss = 1:length(subjIDs)
    if length(stat_file)==1
        stat_file_lh = stat_file{1};
        stat_file_rh = stat_file{1};
    else
        stat_file_lh = stat_file{ss,1};
        stat_file_rh = stat_file{ss,2};
    end
    subjID = subjIDs{ss};

    if use_fsaverage
        surfaceSubjID = 'fsaverage';
    else
        surfaceSubjID = subjID;
    end

    if overlay_path_override
        subj_func_dir = [func_path '/' func_folder '/'];
    else
        subj_func_dir = [func_path subjID '/' func_folder '/'];
    end

    for cc = 1:size(contrast_list,2)
        contrast = contrast_list{ss,cc};
        if overlay_path_override
            lh_overlay = [subj_func_dir '/' lh_analysisName '/' stat_file_lh];
            rh_overlay = [subj_func_dir '/' rh_analysisName '/' stat_file_rh];
        else
            lh_overlay = [subj_func_dir '/' lh_analysisName '/' contrast '/' stat_file_lh];
            rh_overlay = [subj_func_dir '/' rh_analysisName '/' contrast '/' stat_file_rh];
        end

        if ~isempty(lh_label_list) && ~isempty(rh_label_list)
            lh_labelList_curr = lh_label_list{ss,cc};
            rh_labelList_curr = rh_label_list{ss,cc};

            lh_annotList = {};
            rh_annotList = {};
        end

        for h = 1:length(hemis)
            hemi = hemis{h};
            if strcmp(hemi,'lh')
                overlay = lh_overlay;
                if ~isempty(lh_label_list)
                    labelList = lh_labelList_curr;
                    annotList = lh_annotList;
                    label_names = cellfun(@(x) split(x,'.'), lh_labelList_curr, 'UniformOutput', false); % split file name to isolate ROI name
                    label_names = string(cellfun(@(x) x{2}, label_names, 'UniformOutput', false)); % ROI name should be 2nd 
                    inds = cell2mat(cellfun(@(x) find(strcmp(colortable.ROI,x)), label_names, 'UniformOutput', false));
                    labelColors = table2array(colortable(inds, 2:4)); % select ROIs and get last 3 rows for RGB values corresponding to each ROI
                else
                    labelList = {};
                    annotList = {};
                    labelColors = [];
                end
            elseif strcmp(hemi,'rh')
                overlay = rh_overlay;
                if ~isempty(rh_label_list)
                    labelList = rh_labelList_curr;
                    annotList = rh_annotList;
                    label_names = cellfun(@(x) split(x,'.'), rh_labelList_curr, 'UniformOutput', false); % split file name to isolate ROI name
                    label_names = string(cellfun(@(x) x{2}, label_names, 'UniformOutput', false)); % ROI name should be 2nd 
                    inds = cell2mat(cellfun(@(x) find(strcmp(colortable.ROI,x)), label_names, 'UniformOutput', false));
                    labelColors = table2array(colortable(inds, 2:4)); % select ROIs and get last 3 rows for RGB values corresponding to each ROI
                else
                    labelList = {};
                    annotList = {};
                    labelColors = [];
                end
            end
            if no_contrasts
                overlay = '';
            end
            plotFreeviewFunc(subjID,contrast,recon_dir,surfaceSubjID,inflationType,overlay,labelList, ...
                annotList,overlayThreshold,hemi,labelColors,labelOutline,label_opacity,annotOutline,viewAngles,save_dir,ss_suffix, ss_on);
        end
    end
end
disp(['All subjects complete. Run "bash ' save_dir 'run_all_screenshots.sh" from the terminal to start freeview plotting and screenshoting'])

    function plotFreeviewFunc(subjID,contrast,reconDir,surfaceSubjID,inflationType,overlay,labelList, ...
            annotList,overlayThreshold,hemi,labelColors,labelOutline,labelOpacity, annotOutline,viewAngles,saveDir,ss_suffix, ss_on)
        cmd ='freeview';
        cmd = [cmd ' -f ' reconDir '/' surfaceSubjID '/surf/' hemi '.' inflationType ':curvature_method=binary'];

        if ~isempty(overlay)
            cmd = [cmd ':overlay=' overlay ':overlay_threshold=' overlayThreshold];
        end

        for ii = 1:length(labelList)
            label = labelList{ii};
            cmd = [cmd ':label=' label ':label_color=' sprintf('%d,%d,%d', labelColors(ii,:)) ':label_outline=' labelOutline ':label_opacity=' labelOpacity];
        end
        for a = 1:length(annotList)
            annot = annotList{a};
            cmd = [cmd ':annot=' annot ':annot_outline=' annotOutline];
        end

        cmd = [cmd ' -nocursor'];

        % Write the first string to the text file
        fid = fopen([saveDir '/' subjID '_' contrast '_' hemi '.txt'], 'w');
        fprintf(fid, '%s\n', cmd);
        fclose(fid);
        if ss_on
            for v=1:length(viewAngles)
                cmdSS = 'freeview';
                viewAngle=viewAngles{v};
                if strcmp(hemi,'lh')
                    if strcmp(viewAngle,'lateral')
                        cmdSS = [cmdSS];
                    elseif strcmp(viewAngle,'medial')
                        cmdSS = [cmdSS ' -cam Azimuth 180'];
                    elseif strcmp(viewAngle,'posterior')
                        cmdSS = [cmdSS ' -cam Azimuth -90'];
                    end
                elseif strcmp(hemi,'rh')
                    if strcmp(viewAngle,'lateral')
                        cmdSS = [cmdSS ' -cam Azimuth 180'];
                    elseif strcmp(viewAngle,'medial')
                        cmdSS = [cmdSS ' -cam Azimuth 180'];
                    elseif strcmp(viewAngle,'posterior')
                        cmdSS = [cmdSS ' -cam Azimuth 90'];
                    end
                end

                cmdSS = [cmdSS ' -ss ' saveDir '/' subjID '_' contrast '_' viewAngle '_' hemi ss_suffix '.png'];
                fid = fopen([saveDir '/' subjID '_' contrast '_' hemi '.txt'], 'a');
                fprintf(fid, '%s\n', cmdSS);
                fclose(fid);
            end

            fid = fopen([saveDir '/' subjID '_' contrast '_' hemi '.txt'], 'a');
            cmdQuit = 'freeview -quit';
            fprintf(fid, '%s\n', cmdQuit);
            fclose(fid);

        end



        fv_command = ['freeview -cmd ' saveDir '/' subjID '_' contrast '_' hemi '.txt'];
        fid = fopen([saveDir '/run_all_screenshots.sh'], 'a');
        fprintf(fid, '%s\n', fv_command);
        fclose(fid);

        %unix(fv_command);
        %disp(fv_command);
        disp(['Complete ' subjID ' ' contrast]);
    end
end

