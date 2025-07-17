%%
function crop_ppt_fv_images(subjIDs, contrast_list, run_PPT, save_dir, pptTitle)
%FREEVIEW_SCREENSHOTS - Runs through subjs and contrasts
% If matlab is running into errors opening freesurfer, try running this
% with freeview version 6.0 loaded.

% Initialize inputs
assert(nargin >= 2, 'Not enough input arguments (subjIDs, contrast_list, all required)');

if nargin < 3 || isempty(run_PPT)
    run_PPT = true;
end

if nargin < 4 || isempty(save_dir)
    save_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/figures_images/roi_QC_screenshots/';
end

if nargin < 5 || isempty(pptTitle)
    pptTitle='ROI_screenshots.pptx';
end


% Set key variables
hemis = {'lh','rh'};
viewAngles = {'lateral', 'medial'}; % could add dorsal / ventral

no_contrasts = false;
if isempty(contrast_list)
    no_contrasts = true;
    contrast_list = cell(length(subjIDs),1);
    contrast_list(:) = {''};
end

for ss=1:length(subjIDs)
    subjID=subjIDs{ss};
    for cc=1:size(contrast_list,2)
        contrast=contrast_list{ss,cc};
        for h=1:length(hemis)
            hemi=hemis{h};
            for v=1:length(viewAngles)
                viewAngle = viewAngles{v};
                imgFile = [save_dir '/' subjID '_' contrast '_' viewAngle '_' hemi '_fsavg.png']; 
                if strcmp(hemi,'lh')
                    if contains(imgFile,'lateral')
                        cropAmts= [125,125,400,400];
                        %cropAmts= [200,200,550,550]; % top,bottom,left,right
                    elseif contains(imgFile,'medial')
                        cropAmts= [125,125,400,400]; % top,bottom,left,right
                    elseif contains(imgFile,'posterior')
                        cropAmts= [200,200,550,850]; % top,bottom,left,right
                    end
                elseif strcmp(hemi,'rh')
                    if contains(imgFile,'lateral')
                        cropAmts= [125,125,400,400]; % top,bottom,left,right
                    elseif contains(imgFile,'medial')
                        cropAmts= [125,125,400,400]; % top,bottom,left,right
                    elseif contains(imgFile,'posterior')
                        cropAmts= [200,200,850,550]; % top,bottom,left,right
                    end
                end

                cropImageFunc(imgFile,cropAmts(1),cropAmts(2),cropAmts(3),cropAmts(4));
            end

        end
    end
end

if run_PPT
    makePPTFunc(subjIDs,contrast_list,hemis,viewAngles,pptTitle, save_dir);
end


    function cropImageFunc(inputFileName, cropTop, cropBottom, cropLeft, cropRight)
        try
            % Read the image
            imgData = imread(inputFileName);

            % Get the original dimensions
            [origHeight, origWidth, ~] = size(imgData);

            % Calculate the new dimensions after cropping
            newWidth = origWidth - cropLeft - cropRight;
            newHeight = origHeight - cropTop - cropBottom;

            % Check if the cropping values are valid
            if newWidth <= 0 || newHeight <= 0
                error('Cropping values are too large for the image dimensions.');
            end

            % Crop the image
            croppedImg = imgData(cropTop+1:end-cropBottom, cropLeft+1:end-cropRight, :);

            % Save the cropped image to a new file
            newFileName = [inputFileName(1:end-4) '_crop.png'];
            imwrite(croppedImg, newFileName);
        catch
            disp(['Missing ' inputFileName]);
        end
    end

    function makePPTFunc(subjIDs,contrasts,hemis,viewAngles,pptTitle, save_dir)
        import mlreportgen.ppt.*;

        % Create a presentation
        ppt = Presentation(pptTitle);
        open(ppt);

        for ss=1:length(subjIDs)
            subjID = subjIDs{ss};
            for cc=1:size(contrasts,2)
                contrast=contrasts{ss,cc};
                % Add a slide to the presentation
                slide = add(ppt, 'Title and Content');
                replace(slide, 'Title', [subjID ': ' contrast]);

                % Loop through the images and add them to the slide
                for h=1:length(hemis)
                    hemi=hemis{h};
                    for v = 1:length(viewAngles)
                        try
                            viewAngle=viewAngles{v};
                            imgFile = [save_dir '/' subjID '_' contrast '_' viewAngle '_' hemi '_fsavg_crop.png']; 
                            img = Picture(imgFile);

                            % Read the image to get its original dimensions
                            info = imfinfo(imgFile);
                            origWidth = info.Width;
                            origHeight = info.Height;

                            % Set the target width and calculate height to maintain aspect ratio
                            if contains(viewAngle,'posterior')
                                targetWidth = 2.5;
                            else
                                targetWidth = 4; % Target width in inches
                            end
                            targetHeight = (targetWidth / origWidth) * origHeight;

                            % Assign the calculated width and height to the image
                            img.Width = sprintf('%fin', targetWidth);
                            img.Height = sprintf('%fin', targetHeight);

                            % Set the position of the image
                            if strcmp(hemi,'lh')
                                if contains(imgFile,'lateral')
                                    img.X = '0.5in';
                                    img.Y = '2in';
                                elseif contains(imgFile,'medial')
                                    img.X = '0.5in';
                                    img.Y = '4.5in';
                                elseif contains(imgFile,'posterior')
                                    img.X = '4.5in';
                                    img.Y = '2in';
                                end
                            elseif strcmp(hemi,'rh')
                                if contains(imgFile,'lateral')
                                    img.X = '9in';
                                    img.Y = '2in';
                                elseif contains(imgFile,'medial')
                                    img.X = '9in';
                                    img.Y = '4.5in';
                                elseif contains(imgFile,'posterior')
                                    img.X = '6.5in';
                                    img.Y = '2in';
                                end
                            end

                            add(slide, img); % Add image to slide
                        catch
                            disp(['Missing ' imgFile]);
                        end
                    end
                end
            end
        end

        % Close and save the presentation
        close(ppt);

        % Open the presentation using the system's default application
        unix(['libreoffice ' pptTitle]);

    end


end

