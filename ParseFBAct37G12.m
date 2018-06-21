%% Clear out the workspace
clear;
clc;

%% Get the directories and, for each directory, 
% the number of flies, the region that is being imaged, and whether one or two color imaging was performed.
% Go to the root directory
dirRoot = uigetdir('C:\Users\turnerevansd\Downloads\Data','Select the root directory');
cd(dirRoot);
% Specify the directories
numDirs = input('Number of directories? ');
allPathname = cell(numDirs,1);
for dirNow = 1:numDirs
    allPathname{dirNow}.name = uigetdir('C:\Users\turnerevansd\Downloads\Data','Select the directory');
    allPathname{dirNow}.numFlies = input('Number of flies? ');
    allPathname{dirNow}.numColors = input('1 or 2 color imaging? ');
    allPathname{dirNow}.imRegion = input('Region that was imaged? (EB,PB,FB,other) ');
    while ~(strcmp(allPathname{dirNow}.imRegion,'EB') ||...
            strcmp(allPathname{dirNow}.imRegion,'PB') ||...
            strcmp(allPathname{dirNow}.imRegion,'FB') ||...
            strcmp(allPathname{dirNow}.imRegion,'other'))
        allPathname{dirNow}.imRegion = input('Region that was imaged? (EB,PB,other)');
    end
end

%% Average stacks over all R or L turns and compare

% Set a rotational velocity threshold
vRThresh = pi/4;

for dirNow = 1:numDirs
    fileNames = dir(allPathname{dirNow}.name);
    cd(allPathname{dirNow}.name);
    allPathnameNow = strcat(allPathname{dirNow}.name,'\');
    for tifID = 3:length(fileNames)
        tifName = fileNames(tifID).name;
        % Process each trial
        if strcmp(tifName(end-3:end),'.tif')
            % Extract the position and DF/F values
            for moveID = 3:length(fileNames)
                moveName =fileNames(moveID).name;
                tifNameParts = strsplit(tifName,'_');
                moveNameParts = strsplit(moveName,'_');
                % Find the associated position data 
                if (strcmpi(moveName(end-3:end),'.txt') &...
                        strcmp(tifName(end-5:end-4),moveName(end-5:end-4)) &...
                        strcmp(tifName(1:4),moveName(1:4)) &...
                        strcmp(moveNameParts{end-1},tifNameParts{end-1}) & ...
                        strcmp(moveNameParts{end-2},tifNameParts{end-2}))
                    moveName
                    
                    % Load the stacks
                    rotAng = 0;
                    [stackMaxInt, stackMean] = ImDatLoadBigtiff(tifName,allPathnameNow,rotAng);
                    clear stackMean;
                    stackReg = imRegSimple(stackMaxInt, 10); % Image correct the stacks
                    
                    % Load the movement data
                    positionDat = VRDatLoad(moveName,allPathnameNow,0);
                    tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
                    tStep = mean(diff(tSpan));
                    numPlanes = length(positionDat.tFrameGrab)/length(stackMaxInt);
                    minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/numPlanes);
                    maxFG = round(length(positionDat.tFrameGrab)/numPlanes);
                    positionDat.num_planes = numPlanes;
                    positionDat.minFG = minFG;
                    positionDat.maxFG = maxFG;

                    % convert the raw ball position data to position (in case
                    % the visual display is not closed loop with a gain of 1)
                    [posRot, posFor, posLat] = PositionConverter(positionDat);

                    % Match the position data to the framegrab times
                    OffsetRotMatch = MatchData(positionDat.t,positionDat);
                    OffsetRotMatch(:,2) = MatchData(pi/180*positionDat.OffsetRot,positionDat);
                    vR = diff(OffsetRotMatch(:,2))./tStep;
                    OffsetForMatch = MatchData(positionDat.OffsetFor,positionDat);
                    OffsetLatMatch = MatchData(positionDat.OffsetLat,positionDat);
                    PosRotMatch = MatchData(posRot,positionDat);
                    PosForMatch = MatchData(posFor,positionDat);
                    PosLatMatch = MatchData(posLat,positionDat);

                    % Savitsky-Golay filter the rotational velocity
                    vRFilt = sgolayfilt(vR,3,11);
                    
                    % Find the R and L turns
                    RTurn = find(vRFilt > vRThresh);
                    LTurn = find(vRFilt < -vRThresh);
                    
                    % Sort the registered stacks over these times and
                    % average
                    RTurnStack = mean(stackReg(:,:,RTurn),3);
                    LTurnStack = mean(stackReg(:,:,LTurn),3);
                    
                    % Plot the mean stacks and their difference
                    act = figure('units','normalized','outerposition',[0 0 1 1]);

                    % Plot the activity
                    subplot(1,3,1);
                    hold on;
                    imagesc(imrotate(RTurnStack,90));
                    axis equal;
                    axis off;
                    title('R turns');
                    
                    subplot(1,3,2);
                    hold on;
                    imagesc(imrotate(LTurnStack,90));
                    axis equal;
                    axis off;
                    title('L turns');
                    
                    subplot(1,3,3);
                    hold on;
                    imagesc(imrotate(RTurnStack-LTurnStack,90));
                    axis equal;
                    axis off;
                    title('R-L');
                    
                    set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
                    print(act,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\ParseFBAct\DarkCL_',moveNameParts{1},'_',moveNameParts{end-1},'_Trial',moveNameParts{end}(2)),'-dpdf'); 
                    
                    delete(act);

                    save(strcat('C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\ParseFBAct\DarkCL_',moveNameParts{1},'_',moveNameParts{end-1},'_Trial',moveNameParts{end}(2)),'RTurnStack','LTurnStack');
                    
                    clear fullpath stackMaxInt stackReg stackXYfiltMax ROIaveMax positionDat;
                end
            end
        end
    end    
end