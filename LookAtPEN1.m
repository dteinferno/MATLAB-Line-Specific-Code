%% Clear out old data
clear;
clc;

cond = FlyDatLoad(2);

%% Specify the given ROIs
PEN1ROIs = [1:8 11:18];

%% Savitsky-Golay filter all velocities
sgolayOrder = 3;
sgolayFrames = 11;

for condID = 1:length(cond)
    for flyID = 1:cond{condID}.numFlies
        for trialType = 1:3
            if trialType == 1
                trialName = 'Onex';
            elseif trialType == 2
                trialName = 'Onex';
            elseif trialType == 3
                trialName = 'Twox';
            else
                continue;
            end
            for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))
                vFFilt = sgolayfilt(...
                    cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vF,...
                    sgolayOrder,sgolayFrames);
                cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vFFilt = vFFilt;
                vRFilt = sgolayfilt(...
                    cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vRot,...
                    sgolayOrder,sgolayFrames);
                cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vRFilt = vRFilt;
            end
        end
    end
end


%% Plot 2D heading vs. rotational velocity tuning curves for the individual PB ROIs

for condID = 1:length(cond);
    if condID == 1
        imData = 'GROIaveMax';
    else
        imData = 'RROIaveMax';
    end

    % Extract the relevant data
    flyMap = {};
    for flyID=1:cond{condID}.numFlies
        actAll = [];
        vRAll = [];
        OnexPosAll = [];
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.Onex)

            % Extract the appropriate data
            PBDat =  cond{condID}.allFlyData{flyID}.Onex{trialID}.(imData)(PEN1ROIs,1:end-1)-1;
            vR = cond{condID}.allFlyData{flyID}.Onex{trialID}.positionDatMatch.vRFilt;
            OnexPos = cond{condID}.allFlyData{flyID}.Onex{trialID}.positionDatMatch.OffsetRotMatch(1:end-1,2);

            actAll = horzcat(actAll,PBDat);
            vRAll = vertcat(vRAll,vR);
            OnexPosAll = vertcat(OnexPosAll,OnexPos);
        end
        flyMap{flyID}.actAll = actAll;
        flyMap{flyID}.vRAll = vRAll;
        flyMap{flyID}.OnexPosAll = OnexPosAll;
    end

    % Create a 2D matrix and tuning curves out of the data
    headingEdges = linspace(-pi,pi,17);
    headingCents = headingEdges(1:end-1);
    headingCents = headingCents+0.5*mean(diff(headingCents));
    vREdges = linspace(-3*pi,3*pi,49);
    vRCents = vREdges(1:end-1);
    vRCents = vRCents+0.5*mean(diff(vRCents));

    for flyID=1:cond{condID}.numFlies
        % Initialize the objects and arrays
        mapMatrixAll = {};
        mapMatrixMean = zeros(16,length(headingEdges)-1,length(vREdges)-1);
        headingTuneAll = {};
        headingTuneMean = zeros(16,length(headingEdges)-1);
        vRTuneAll = {};
        vRTuneMean = zeros(16,length(vREdges)-1);

        % Bin the heading and velocity data
        headBin = discretize(flyMap{flyID}.OnexPosAll,headingEdges);
        vRBin = discretize(flyMap{flyID}.vRAll,vREdges);
        minN = 10;

        % Step through the ROIs
        for ROINow = 1:16
            % Populate the 2D array
            for headingStep = 1:length(headingEdges)-1
                for vRStep = 1:length(vREdges)-1
                    mapMatrixAll{ROINow}.heading{headingStep}.vR{vRStep} = [];
                end
            end
            for tPt = 1:length(flyMap{flyID}.vRAll)
                mapMatrixAll{ROINow}.heading{headBin(tPt)}.vR{vRBin(tPt)} = ...
                    vertcat(mapMatrixAll{ROINow}.heading{headBin(tPt)}.vR{vRBin(tPt)},...
                    flyMap{flyID}.actAll(ROINow,tPt));
            end
            for headingStep = 1:length(headingEdges)-1
                for vRStep = 1:length(vREdges)-1
                    mapMatrixMean(ROINow,headingStep,vRStep) = ...
                        mean(mapMatrixAll{ROINow}.heading{headingStep}.vR{vRStep});
                end
            end

            % Populate the linear arrays
            for headingStep = 1:length(headingEdges)-1
                headingTuneAll{ROINow}.heading{headingStep} = [];
            end
            for vRStep = 1:length(vREdges)-1
                vRTuneAll{ROINow}.vR{vRStep} = [];
            end
            for tPt = 1:length(flyMap{flyID}.vRAll)
                headingTuneAll{ROINow}.heading{headBin(tPt)} = ...
                    vertcat(headingTuneAll{ROINow}.heading{headBin(tPt)},...
                    flyMap{flyID}.actAll(ROINow,tPt));
                vRTuneAll{ROINow}.vR{vRBin(tPt)} = ...
                    vertcat(vRTuneAll{ROINow}.vR{vRBin(tPt)},...
                    flyMap{flyID}.actAll(ROINow,tPt));
            end
            for headingStep = 1:length(headingEdges)-1
                if length(headingTuneAll{ROINow}.heading{headingStep}) > minN
                    headingTuneMean(ROINow,headingStep) =...
                         mean(headingTuneAll{ROINow}.heading{headingStep});
                else
                    headingTuneMean(ROINow,headingStep) = NaN;
                end
            end
            for vRStep = 1:length(vREdges)-1
                if length(vRTuneAll{ROINow}.vR{vRStep}) > minN
                    vRTuneMean(ROINow,vRStep) =...
                         mean(vRTuneAll{ROINow}.vR{vRStep});
                else
                    vRTuneMean(ROINow,vRStep) = NaN;
                end
            end
        end
        flyMap{flyID}.mapMatrixAll = mapMatrixAll;
        flyMap{flyID}.mapMatrixMean = mapMatrixMean;
        flyMap{flyID}.headingTuneAll = headingTuneAll;
        flyMap{flyID}.headingTuneMean = headingTuneMean;
        flyMap{flyID}.vRTuneAll = vRTuneAll;
        flyMap{flyID}.vRTuneMean = vRTuneMean;
    end

    for flyID=1:cond{condID}.numFlies
        vRHeadingMapL = figure('units','normalized','outerposition',[0 0 1 1]);
        for ROINow = 1:8
            if ROINow <=4
                subplot(6,12,[1:2 13:14]+3*(ROINow-1));
                imagesc(vRCents,headingCents,flipud(squeeze(flyMap{flyID}.mapMatrixMean(ROINow,:,:))));
                axis off;
                title(strcat('ROI #',num2str(ROINow)));
                if ROINow == 1
                    text(-0.25,-pi-0.5,strcat('fly #',num2str(flyID),'-L PB'),'FontSize',14);
                end

                subplot(6,12,[25:26]+3*(ROINow-1));
                plot(vRCents,flyMap{flyID}.vRTuneMean(ROINow,:),'Color',[1 0.5 0]);
                xlabel('vR (rad/s)');
                xlim([vRCents(1) vRCents(end)]);
                ylabel('DF/F');

                subplot(6,12,[3 15]+3*(ROINow-1));
                plot(flyMap{flyID}.headingTuneMean(ROINow,:),headingCents,'Color',[0 0 1]);
                xlabel('DF/F');            
                ylabel('heading (rad)');
                ylim([headingCents(1) headingCents(end)]);

            else
                subplot(6,12,[37:38 49:50]+3*(ROINow-5))
                imagesc(vRCents,headingCents,flipud(squeeze(flyMap{flyID}.mapMatrixMean(ROINow,:,:))));
                axis off;
                title(strcat('ROI #',num2str(ROINow)));

                subplot(6,12,[61:62]+3*(ROINow-5));
                plot(vRCents,flyMap{flyID}.vRTuneMean(ROINow,:),'Color',[1 0.5 0]);
                xlabel('vR (rad/s)');
                xlim([vRCents(1) vRCents(end)]);
                ylabel('DF/F');

                subplot(6,12,[39 51]+3*(ROINow-5));
                plot(flyMap{flyID}.headingTuneMean(ROINow,:),headingCents,'Color',[0 0 1]);
                xlabel('DF/F');            
                ylabel('heading (rad)');
                ylim([headingCents(1) headingCents(end)]);
            end
        end
        set(vRHeadingMapL,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
        print(vRHeadingMapL,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\TurningPaper\PEN1-PB\Figures\vRHeadingMap_LPB_',cond{condID}.allFlyData{flyID}.ID),'-dpdf');
        delete(vRHeadingMapL);

        vRHeadingMapR = figure('units','normalized','outerposition',[0 0 1 1]);
        for ROINow = 1:8
            if ROINow <=4
                subplot(6,12,[1:2 13:14]+3*(ROINow-1));
                imagesc(vRCents,headingCents,flipud(squeeze(flyMap{flyID}.mapMatrixMean(ROINow+4,:,:))));
                axis off;
                title(strcat('ROI #',num2str(ROINow+4)));
                if ROINow == 1
                    text(-0.25,-pi-0.5,strcat('fly #',num2str(flyID),'-R PB'),'FontSize',14);
                end

                subplot(6,12,[25:26]+3*(ROINow-1));
                plot(vRCents,flyMap{flyID}.vRTuneMean(ROINow+4,:),'Color',[1 0.5 0]);
                xlabel('vR (rad/s)');
                xlim([vRCents(1) vRCents(end)]);
                ylabel('DF/F');

                subplot(6,12,[3 15]+3*(ROINow-1));
                plot(flyMap{flyID}.headingTuneMean(ROINow+4,:),headingCents,'Color',[0 0 1]);
                xlabel('DF/F');            
                ylabel('heading (rad)');
                ylim([headingCents(1) headingCents(end)]);

            else
                subplot(6,12,[37:38 49:50]+3*(ROINow-5))
                imagesc(vRCents,headingCents,flipud(squeeze(flyMap{flyID}.mapMatrixMean(ROINow+4,:,:))));
                axis off;
                title(strcat('ROI #',num2str(ROINow+4)));

                subplot(6,12,[61:62]+3*(ROINow-5));
                plot(vRCents,flyMap{flyID}.vRTuneMean(ROINow+4,:),'Color',[1 0.5 0]);
                xlabel('vR (rad/s)');
                xlim([vRCents(1) vRCents(end)]);
                ylabel('DF/F');

                subplot(6,12,[39 51]+3*(ROINow-5));
                plot(flyMap{flyID}.headingTuneMean(ROINow+4,:),headingCents,'Color',[0 0 1]);
                xlabel('DF/F');            
                ylabel('heading (rad)');
                ylim([headingCents(1) headingCents(end)]);
            end
        end
        set(vRHeadingMapR,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
        print(vRHeadingMapR,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\TurningPaper\PEN1-PB\Figures\vRHeadingMap_RPB_',cond{condID}.allFlyData{flyID}.ID),'-dpdf');
        delete(vRHeadingMapR);

    end
end
