function cond = FlyDatLoad(numColors, region)
%%returns cond, an object storing some data for each condition
    %% Prompt the user for the root directory
    %dirRoot = uigetdir('C:/Users/turnerevansd/Documents','Select the root directory');
    dirRoot = '/Users/loaner/Documents/imaging/Data_Dan/shi'; %Kris

    %% Prompt the user for the number of conditions
    numConditions = input('# of conditions?');
    %numConditions = 1

    %% Ask the user for the name of each condition and the directories
    % associated with that condition
    dirL = 0;
    cond = cell(numConditions,1); %new cell for each condition
    for condStep = 1:numConditions
        cond{condStep}.name = input(strcat('name of cond #',num2str(condStep),'?')); %add name
        %cond{condStep}.name = region %Kris
        dirL = input(strcat('Total # of dirs for cond:',cond{condStep}.name,'?')); %add dir count
        %dirL = 1 %Kris
        cond{condStep}.dirs = cell(dirL,1); %add dirs
        for dirStep = 1:dirL
            cond{condStep}.dirs{dirStep} = uigetdir(dirRoot,strcat('Specify dir # ',num2str(dirStep)));
            %cond{condStep}.dirs{dirStep}='/Users/loaner/Documents/Imaging/Data_Dan/EPG_G_GE_R_1' %Kris
        end
    end

    %% Find the name for each fly and load the data
    % step through the conditions
    for condStep =1:numConditions
        allFlyData{1}.ID = 'Empty'; % a cell to hold the data for each fly
        numFlies = 1;

        % step through the directories to load the pre-processed data for each
        % fly
        for dirStep = 1:length(cond{condStep}.dirs) %for each directory
            %Can e.g. have multiple directories with the same kind of data
            %but from different dates
            fileNames = dir(cond{condStep}.dirs{dirStep});
            dirParts = strsplit(cond{condStep}.dirs{dirStep},'/');
            for fileID = 3:length(fileNames)
                fileName = fileNames(fileID).name;
                % check to see if the file is the pre-processed data
                if strcmpi(fileName(end-2:end),'mat') & ~strcmpi(fileName(end-7:end-4),'ROIs')
                    nameParts = strsplit(fileName,'_');
                    % specify the fly name
%                     flyName = strcat(dirParts{end},'-',nameParts{1},'-',nameParts{4});
                    flyName = strcat(dirParts{end},'-',nameParts{1},'-',nameParts{3});
                    for flyStep = 1:length(allFlyData)
                        if strcmp(flyName,allFlyData{flyStep}.ID)
                            break;
                        end
                        if flyStep == length(allFlyData)
                            allFlyData{numFlies}.ID = flyName;
                            numFlies = numFlies + 1;
                        end
                    end

                    % Get the position data, ROI data, and trial ID for the
                    % given fly and trial
                    load(strcat(cond{condStep}.dirs{dirStep},'/',fileName));
                    strcat(cond{condStep}.dirs{dirStep},'/',fileName)
                    
                    if numColors == 1
                        ROIaveMaxFilt = zeros(size(ROIaveMax));
                        if ndims(ROIaveMaxFilt) == 2
                            ROIaveMaxFilt = ROIaveMax;
                        else
                            for ROIgp=1:size(ROIaveMaxFilt,1)
                                ROIaveMaxFilt(ROIgp,:,:) = squeeze(ROIaveMax(ROIgp,:,:));
                            end
                        end
                    else
                        RROIaveMaxFilt = zeros(size(RROIaveMax));
                        GROIaveMaxFilt = zeros(size(GROIaveMax));
                        if ndims(RROIaveMaxFilt) == 2
                            RROIaveMaxFilt = RROIaveMax; %have data per time per G
                        else
                            for ROIgp=1:size(RROIaveMaxFilt,1)
                                RROIaveMaxFilt(ROIgp,:,:) = squeeze(RROIaveMax(ROIgp,:,:));
                            end
                        end
                        if ndims(GROIaveMaxFilt) == 2
                            GROIaveMaxFilt = GROIaveMax;
                        else
                            for ROIgp=1:size(GROIaveMaxFilt,1)
                                GROIaveMaxFilt(ROIgp,:,:) = squeeze(GROIaveMax(ROIgp,:,:));
                            end
                        end
                    end

                    % Match the behavior to the imaging
                    % First, find the experimental parameters
                    if numColors == 1
                        num_planes = round(length(positionDat.tFrameGrab)/length(ROIaveMax));
                    else
                        num_planes = round(length(positionDat.tFrameGrab)/length(RROIaveMax));
                    end
                    tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
                    minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
                    maxFG = round(length(positionDat.tFrameGrab)/num_planes);
                    positionDat.num_planes = num_planes;
                    positionDat.minFG = minFG;
                    positionDat.maxFG = maxFG;

                    % convert the raw ball position data to position (in case
                    % the visual display is not closed loop with a gain of 1)
                    [posRot, posFor, posLat] = PositionConverter(positionDat);

                    %all of this is similar to the other script
                    
                    % Match the position data to the framegrab times
                    OffsetRotMatch = MatchData(positionDat.t,positionDat);
                    OffsetRotMatch(:,2) = MatchData(pi/180*positionDat.OffsetRot,positionDat);
                    OffsetForMatch = MatchData(positionDat.OffsetFor,positionDat);
                    OffsetLatMatch = MatchData(positionDat.OffsetLat,positionDat);
                    PosRotMatch = MatchData(posRot,positionDat);
                    PosForMatch = MatchData(posFor,positionDat);
                    PosLatMatch = MatchData(posLat,positionDat);
                    Closed = MatchData(positionDat.closed,positionDat);
                    Direction = MatchData(positionDat.direction,positionDat);
                    Trans = MatchData(positionDat.trans,positionDat);
                    if isfield(positionDat,'olgain')
                        OLGain = MatchData(positionDat.olgain,positionDat);
                    end
                    if isfield(positionDat,'clgain')
                        CLGain = MatchData(positionDat.clgain,positionDat);
                    end
                    
%                     if isfield(positionDat,'tStim')
%                         tStim = zeros(positionDat.maxFG-positionDat.minFG+1,1);
%                         for stim = 1:length(positionDat.tStim)
%                             tMatch = find(positionDat.tFrameGrab > positionDat.tStim(stim));
%                             volNow = ceil(tMatch(1)./positionDat.num_planes) - ...
%                                 positionDat.minFG+1;
%                             tStim(volNow) = 1;
%                         end
%                     end
                        

                    % Unwrap the rotation and calculate the velocities
                    OffsetRotUnwrap = UnWrap(pi/180*PosRotMatch,2,0);
                    vRot = diff(OffsetRotUnwrap)./mean(diff(OffsetRotMatch(:,1)));
                    vF = sqrt(diff(PosForMatch).^2+diff(PosLatMatch).^2)./mean(diff(OffsetRotMatch(:,1)));

                    % populate the positionDatMatch cell
                    positionDatMatch = {};
                    positionDatMatch.OffsetRotMatch = OffsetRotMatch;
                    positionDatMatch.OffsetForMatch = OffsetForMatch;
                    positionDatMatch.OffsetLatMatch = OffsetLatMatch;
                    positionDatMatch.vRot = vRot;
                    positionDatMatch.vF = vF;
                    positionDatMatch.minFG = minFG;
                    positionDatMatch.maxFG = maxFG;
                    positionDatMatch.PosRotMatch = PosRotMatch;
                    positionDatMatch.PosForMatch = PosForMatch;
                    positionDatMatch.PosLatMatch = PosLatMatch;
                    positionDatMatch.Closed = Closed;
                    positionDatMatch.Direction = Direction;
                    positionDatMatch.Trans = Trans;
                    if isfield(positionDat,'olgain')
                        positionDatMatch.OLGain = OLGain;
                    end
                    if isfield(positionDat,'clgain')
                        positionDatMatch.CLGain = CLGain;
                    end
%                     if isfield(positionDat,'tStim')
%                         positionDatMatch.tStim = tStim;
%                     end

                    % load the matched position data and imaging data into the
                    % allFlyData cell
                    for flyID = 1:length(allFlyData)
                        if strcmp(flyName,allFlyData{flyID}.ID)
                            fileParts = strsplit(fileName,'_');
                            trialType = fileParts{end-1};
                            if strcmp(trialType(1),'1')
                                trialType = strcat('One',trialType(2:end));
                            elseif strcmp(trialType(1),'2')
                                trialType = strcat('Two',trialType(2:end));
                            end
                            if ismember('30C',fileParts)
                                trialType = strcat(trialType,'_30C'); %high T; shibire. NO endocytosis
                            end
                            trialID = str2num(fileName(end-5:end-4));

                            allFlyData{flyID}.(trialType){trialID}.positionDatMatch = positionDatMatch;
                            trialType
                            trialID
                            if numColors == 1
                                if ndims(ROIaveMaxFilt) == 2
                                    allFlyData{flyID}.(trialType){trialID}.ROIaveMax = ROIaveMaxFilt(:,minFG:maxFG);
                                else
                                    allFlyData{flyID}.(trialType){trialID}.ROIaveMax = ROIaveMaxFilt(:,:,minFG:maxFG);                                    
                                end
                            else
                                if ndims(RROIaveMaxFilt) == 2
                                    allFlyData{flyID}.(trialType){trialID}.RROIaveMax = RROIaveMaxFilt(:,minFG:maxFG);
                                else
                                    allFlyData{flyID}.(trialType){trialID}.RROIaveMax = RROIaveMaxFilt(:,:,minFG:maxFG);
                                end
                                if ndims(GROIaveMaxFilt) == 2
                                    allFlyData{flyID}.(trialType){trialID}.GROIaveMax = GROIaveMaxFilt(:,minFG:maxFG);
                                else
                                    allFlyData{flyID}.(trialType){trialID}.GROIaveMax = GROIaveMaxFilt(:,:,minFG:maxFG);
                                end
                            end
                        end
                    end
                end
            end
        end

        cond{condStep}.allFlyData = allFlyData;
        cond{condStep}.numFlies = numFlies-1;
        clear allFlyData;
    end

end