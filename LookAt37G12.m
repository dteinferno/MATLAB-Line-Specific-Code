%% Clear out old data
clear;
clc;

cond = FlyDatLoad(1);

%% Set velocity thresholds
vFThresh = 0.1;
vRThresh = pi/32;

%% Savitsky-Golay filter all velocities
sgolayOrder = 3;
sgolayFrames = 11;

for condID = 1:length(cond)
    for flyID = 1:cond{condID}.numFlies
        for trialType = 1:2
            if condID == 1 & trialType == 1
                trialName = 'Dark';
            elseif condID == 1 & trialType == 2
                trialName = 'Stripe';
            elseif condID == 2 & trialType == 1
                trialName = 'Flow';
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

%% Plot an example for the stripe case

flyEx = 7;
trialEx = 1;

if 0
    % Plot the stack mean
    imagePathname = 'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\DarkCL\20180328\';
    imageFilename = 'Fly2_5day_6fx37G12_Stripe_00001.tif';
    rotAng = -70;
    [stackMean, stackMaxInt] = ImDatLoadBigtiff(imageFilename,imagePathname,rotAng);
    meanIm = (mean(stackMaxInt,3));
    meanIm(find(meanIm == 0)) = min(min(meanIm));
    meanIm = (meanIm - min(min(meanIm)))./(max(max(meanIm))-min(min(meanIm)));
    meanIm(find(meanIm < 0.3)) = 0.3;
    meanIm = (meanIm - min(min(meanIm)))./(max(max(meanIm))-min(min(meanIm)));
end

exAct = figure('units','normalized','outerposition',[0 0 1 1]);

pltIm = zeros(size(meanIm,1),size(meanIm,2),3);
pltIm(:,:,2) = meanIm;
subplot(3,3,1);
hold on;
image(flipud(fliplr(pltIm)));

% Plot the ROIs
load('C:\Users\turnerevansd\Documents\RawAnalysis\37G12\DarkCL\20180328\Fly2_5day_6fx37G12_ReferenceROIs');
for ROINow = 1:length(position)
    coordsNow = position{ROINow};
    coordsNow(:,1) = size(pltIm,1) - coordsNow(:,1);
    coordsNow(:,2) = size(pltIm,2) - coordsNow(:,2);
    for lineNow = 1:size(coordsNow,1);
        if lineNow < size(coordsNow,1);
            line([coordsNow(lineNow,1) coordsNow(lineNow+1,1)],...
                [coordsNow(lineNow,2) coordsNow(lineNow+1,2)],...
                'Color','w','LineStyle','--');
        else
            line([coordsNow(lineNow,1) coordsNow(1,1)],...
                [coordsNow(lineNow,2) coordsNow(1,2)],...
                'Color','w','LineStyle','--');
        end
    end
end

% Restrict the range
xlim([30 280]);
ylim([75 230]);
axis equal;
axis off;

% Find the dark periods
DarkPer = find(cond{1}.allFlyData{flyEx}.Stripe{trialEx}.positionDatMatch.Direction == 0);

% Find the closed loop period
CLPer = find(cond{1}.allFlyData{flyEx}.Stripe{trialEx}.positionDatMatch.Closed > 0);
CLPer = setdiff(CLPer,DarkPer);

% Plot the FB activity over time
tPts = cond{1}.allFlyData{flyEx}.Stripe{trialEx}.positionDatMatch.OffsetRotMatch(CLPer,1);
tPts = tPts - tPts(1);
ROIs = [1:8];
FBAct = cond{1}.allFlyData{flyEx}.Stripe{trialEx}.ROIaveMax(:,CLPer)-1;

subplot(3,3,[2:3]);
imagesc(tPts,ROIs,FBAct);
ylabel('ROI #');
title('FB activity');
colormap(brewermap(64, 'Blues'));
colorbar;
caxis([0 1.5]);

% Plot the max activity vs. the forward velocity
vF = cond{1}.allFlyData{flyEx}.Stripe{trialEx}.positionDatMatch.vFFilt(CLPer(1:end-1));
maxF = max(FBAct);
subplot(3,3,[5:6]);
hold on;
plot(tPts(1:end-1),vF,'Color',[1 0.75 0]);
plot(tPts(1:end-1),maxF(2:end),'Color','g');
legend({'vF','max DF/F'});
legend('boxoff');
xlim([tPts(1) tPts(end)]);
ylim([0 1.5]);
ylabel('vF (cm/s), DF/F');
colorbar;

% Plot the heading vs. the stripe position
% Find the PVA
num_ROIs = size(FBAct,1);
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';
clear meanAngRaw;
clear meanIntRaw;
for ts = 1:size(FBAct,2)
    meanAngRaw(ts) = circ_mean(angsraw,...
        squeeze(FBAct(:,ts)));
    meanIntRaw(ts) = circ_r(angsraw,...
        squeeze(FBAct(:,ts)));
end

% Find the stripe position
stripePos = cond{1}.allFlyData{flyEx}.Stripe{trialEx}.positionDatMatch.OffsetRotMatch(CLPer,2);

% Remove jumps so that I can plot properly
PVAPlt = meanAngRaw;
stripePosPlt = stripePos;
jumpPVA = [];
jumpStripe = [];
for jump = 1:length(PVAPlt)-1
    if abs(PVAPlt(jump+1)-PVAPlt(jump)) > pi/2
        jumpPVA = vertcat(jumpPVA,jump+1);
    end
    if abs(stripePosPlt(jump+1)-stripePosPlt(jump)) > pi/2
        jumpStripe = vertcat(jumpStripe,jump+1);
    end
end
PVAPlt(jumpPVA) = NaN;
stripePosPlt(jumpStripe) = NaN;


subplot(3,3,[8:9]);
hold on;
plot(tPts,stripePosPlt,'b');
plot(tPts,PVAPlt,'k');
legend({'heading','PVA'});
legend('boxoff');
xlim([tPts(1) tPts(end)]);
ylim([-pi pi]);
xlabel('time (s)');
ylabel('position (rad)');
colorbar;

% Find a fit between the DF max and the velocity
subplot(3,3,4);
hold on;
flyMov = find(vF>vFThresh);
scatter(vF(flyMov),maxF(flyMov+1),20,'g','filled');
alpha(0.2);
[p,S] = polyfit(vF(flyMov),maxF(flyMov+1)',1);
vFRng = linspace(0,0.5,21);
DFFit = polyval(p,vFRng);
plot(vFRng,DFFit,'g');
xlabel('vF (cm/s)');
ylabel('DF/F');

% Plot the real and derived position
XReal = cond{1}.allFlyData{flyEx}.Stripe{trialEx}.positionDatMatch.PosForMatch;
YReal = cond{1}.allFlyData{flyEx}.Stripe{trialEx}.positionDatMatch.PosLatMatch;
pReal = polyfit(XReal,YReal,1);
XInner = zeros(length(vF),1);
YInner = zeros(length(vF),1);
tStep = mean(diff(tPts));
for posStep = 2:length(XInner);
    XInner(posStep) = XInner(posStep-1)+tStep*cos(meanAngRaw(posStep-1))*(maxF(posStep-1)-p(2))./p(1);
    YInner(posStep) = YInner(posStep-1)+tStep*sin(meanAngRaw(posStep-1))*(maxF(posStep-1)-p(2))./p(1);
end
pInner = polyfit(XInner,YInner,1);
angDiff = atan(pReal(1))-atan(pInner(1));
InnerRot = [cos(angDiff), -sin(angDiff);sin(angDiff), cos(angDiff)]*[XInner YInner]';

subplot(3,3,7);
hold on;
plot(XReal,YReal,'b');
plot(InnerRot(1,:),InnerRot(2,:),'r');
legend({'real position','activity derived position'});
legend('boxoff');
xlabel('x pos (cm)');
ylabel('y pos (cm)');
axis equal;

set(exAct,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(exAct,'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\Example','-dpdf');

%% Plot autocorrelation functions for the dark and stripe cases for the filtered velocities

% Set the number of point around which to consider the autocorrelation
numPts = 60;

% Calculate the autocorrelations
% Step through the flies
for flyID = 1:cond{1}.numFlies
    
    % Step through the trial types
    for trialType = 1:2
        if trialType == 1
            trialName = 'Dark';
        else
            trialName = 'Stripe';
        end
        for trialID=1:length(cond{1}.allFlyData{flyID}.(trialName))

            tStep = mean(diff(cond{1}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1)));
            autoCvF = zeros(numPts,1);
            autoCvR = zeros(numPts,1);
            for lag=1:numPts
                vF = cond{1}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vFFilt;
                vR = cond{1}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vRFilt;
                maxVals = max(cond{1}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax);
                lagCCsF = corrcoef(...
                    vF(floor(numPts/2)-1:end-ceil(numPts/2)-1),...
                    maxVals(numPts-lag+1:end-lag));
                lagCCsR = corrcoef(...
                    vR(floor(numPts/2)-1:end-ceil(numPts/2)-1),...
                    maxVals(numPts-lag+1:end-lag));
                autoCvF(lag) = lagCCsF(2,1);
                autoCvR(lag) = lagCCsR(2,1);
            end
            cond{1}.allFlyData{flyID}.(trialName){trialID}.autoCvF = autoCvF;
            cond{1}.allFlyData{flyID}.(trialName){trialID}.autoCvR = autoCvR;
        end
    end
end
 %% Plot the autocorrelations - per trial
autoCPltPerTrial = figure('units','normalized','outerposition',[0 0 1 1]);

% Step throught the flies
for flyID = 1:cond{1}.numFlies
    
    % Step through the trial types
    for trialType = 1:2
        if trialType == 1
            trialName = 'Dark';
            pltColor = 'k';
        else
            trialName = 'Stripe';
            pltColor = 'b';
        end
        for trialID=1:length(cond{1}.allFlyData{flyID}.(trialName))

            tStep = mean(diff(cond{1}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1)));
            tPts = tStep*([1:numPts]-31);
            autoCvF = cond{1}.allFlyData{flyID}.(trialName){trialID}.autoCvF;
            autoCvR = cond{1}.allFlyData{flyID}.(trialName){trialID}.autoCvR;
            
            subplot(2,cond{1}.numFlies,flyID);
            hold on;
            plot(tPts,autoCvF,pltColor);
            xlim([tPts(1) tPts(end)]);
            ylim([-0.4 1]);
            line([tPts(1) tPts(end)],[0 0],'Color','k');
            line([0 0],[-0.4 1],'Color','k');
            if flyID == 1 & trialID == 1
                ylabel('correlation');
                text(tPts(1)-2,1.1,'vF','FontSize',14);
                line([tPts(end)-3 tPts(end)-2],[0.8 0.8],'Color','k');
                text(tPts(end)-1.5,0.8,'dark trials');
                line([tPts(end)-3 tPts(end)-2],[0.6 0.6],'Color','b');
                text(tPts(end)-1.5,0.6,'CL stripe trials');
            end
            
            subplot(2,cond{1}.numFlies,cond{1}.numFlies+flyID);
            hold on;
            plot(tPts,autoCvR,pltColor);
            xlim([tPts(1) tPts(end)]);
            ylim([-0.4 1]);
            line([tPts(1) tPts(end)],[0 0],'Color','k');
            line([0 0],[-0.4 1],'Color','k');
            xlabel('time (s)');
            if flyID == 1 & trialID == 1
                ylabel('correlation');
                text(tPts(1)-2,1.1,'vR','FontSize',14);
            end
        end
    end
end
 %% Plot the autocorrelations - averaged over all trials
autoCPltAve = figure('units','normalized','outerposition',[0 0 1 1]);

% Step throught the flies
for flyID = 1:cond{1}.numFlies

    clear allCvF allCvR meanCvF meanCvR stdCvR stdCvF;
    % Step through the trial types
    for trialType = 1:2
        if trialType == 1
            trialName = 'Dark';
            pltColor = 'k';
        else
            trialName = 'Stripe';
            pltColor = 'b';
        end
        allCvF{trialType} =  [];
        allCvR{trialType} =  [];
        for trialID = 1:length(cond{1}.allFlyData{flyID}.(trialName))

            autoCvF = cond{1}.allFlyData{flyID}.(trialName){trialID}.autoCvF;
            allCvF{trialType} = horzcat(allCvF{trialType},autoCvF);
            autoCvR = cond{1}.allFlyData{flyID}.(trialName){trialID}.autoCvR;
            allCvR{trialType} = horzcat(allCvR{trialType},autoCvR);
            
        end
        tStep = mean(diff(cond{1}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1)));
        tPts = tStep*([1:numPts]-31);
    end
    
    % Calculate the mean and SD for the trails
    meanCvF{1} = mean(allCvF{1},2)';
    meanCvF{2} = mean(allCvF{2},2)';
    stdCvF{1} = std(allCvF{1},[],2)';
    stdCvF{2} = std(allCvF{1},[],2)';
    cond{1}.allFlyData{flyID}.meanCvF = meanCvF;
    cond{1}.allFlyData{flyID}.stdCvF = stdCvF;
    
    subplot(2,cond{1}.numFlies,flyID);
    hold on;
    plot(tPts,meanCvF{1},'k');
    p1 = patch('XData',[tPts fliplr(tPts)],...
        'YData',[(meanCvF{1}-stdCvF{1}) fliplr(meanCvF{1}+stdCvF{1})],...
        'FaceColor','k','FaceAlpha',0.4);
    uistack(p1,'bottom');
    plot(tPts,meanCvF{2},'b');
    p2 = patch('XData',[tPts fliplr(tPts)],...
        'YData',[(meanCvF{2}-stdCvF{2}) fliplr(meanCvF{2}+stdCvF{2})],...
        'FaceColor','b','FaceAlpha',0.4);
    uistack(p2,'bottom');
    
    xlim([tPts(1) tPts(end)]);
    ylim([-0.4 1]);
    line([tPts(1) tPts(end)],[0 0],'Color','k');
    line([0 0],[-0.4 1],'Color','k');
    if flyID == 1
        ylabel('correlation');
        text(tPts(1)-2,1.1,'vF','FontSize',14);
        line([tPts(end)-3 tPts(end)-2],[0.8 0.8],'Color','k');
        text(tPts(end)-1.5,0.8,'dark mean');
        line([tPts(end)-3 tPts(end)-2],[0.6 0.6],'Color','b');
        text(tPts(end)-1.5,0.6,'CL stripe mean');
    end

    % Calculate the mean and SD for the trials
    meanCvR{1} = mean(allCvR{1},2)';
    meanCvR{2} = mean(allCvR{2},2)';
    stdCvR{1} = std(allCvR{1},[],2)';
    stdCvR{2} = std(allCvR{1},[],2)';
    cond{1}.allFlyData{flyID}.meanCvR = meanCvR;
    cond{1}.allFlyData{flyID}.stdCvR = stdCvR;
    
    subplot(2,cond{1}.numFlies,cond{1}.numFlies+flyID);
    hold on;
    plot(tPts,meanCvR{1},'k');
    p1 = patch('XData',[tPts fliplr(tPts)],...
        'YData',[(meanCvR{1}-stdCvR{1}) fliplr(meanCvR{1}+stdCvR{1})],...
        'FaceColor','k','FaceAlpha',0.4);
    uistack(p1,'bottom');
    plot(tPts,meanCvR{2},'b');
    p2 = patch('XData',[tPts fliplr(tPts)],...
        'YData',[(meanCvR{2}-stdCvR{2}) fliplr(meanCvR{2}+stdCvR{2})],...
        'FaceColor','k','FaceAlpha',0.4);
    uistack(p2,'bottom');
    xlim([tPts(1) tPts(end)]);
    ylim([-0.4 1]);
    line([tPts(1) tPts(end)],[0 0],'Color','k');
    line([0 0],[-0.4 1],'Color','k');
    xlabel('time (s)');
    if flyID == 1
        ylabel('vR correlation');
        text(tPts(1)-2,1.1,'vR','FontSize',14);
    end
end
 %% Plot the peak correlation coefficients
autoCPltPeak = figure('units','normalized','outerposition',[0 0 1 1]);

% Pick an example fly to showcase
flyEx = 1;

% Step throught the flies
for flyID = 1:cond{1}.numFlies

    meanCvF = cond{1}.allFlyData{flyID}.meanCvF;
    meanCvR = cond{1}.allFlyData{flyID}.meanCvR;
    
    subplot(2,2,[2 4]);
    hold on;
    scatter(1,max(meanCvF{1}),200,[1 0.75 0.5],'filled');
    alpha(0.6);
    scatter(2,max(meanCvR{1}),200,[0 0.5 1],'filled');
    alpha(0.6);
    line([1 2],[max(meanCvF{1}) max(meanCvR{1})],'Color',[0.9 0.9 0.9]);
    if flyID == flyEx
        scatter(1,max(meanCvF{1}),200,'k');
        scatter(2,max(meanCvR{1}),200,'k');
    end
        
    
    xlim([0 3]);
    ylim([0 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'vF', 'vR'});
    ylabel('correlation');
    
    if flyID == flyEx
        subplot(2,2,1);
        hold on;
        plot(tPts,meanCvF{1},'Color',[1 0.75 0.5]);
        p1 = patch('XData',[tPts fliplr(tPts)],...
            'YData',[(meanCvF{1}-stdCvF{1}) fliplr(meanCvF{1}+stdCvF{1})],...
            'FaceColor',[1 0.75 0.5],'FaceAlpha',0.4);
        uistack(p1,'bottom');

        xlim([tPts(1) tPts(end)]);
        ylim([-0.4 1]);
        line([tPts(1) tPts(end)],[0 0],'Color','k');
        line([0 0],[-0.4 1],'Color','k');
        ylabel('correlation');
        title('vF vs. DF/F correlation');

        subplot(2,2,3);
        hold on;
        plot(tPts,meanCvR{1},'Color',[0 0.5 1]);
        p1 = patch('XData',[tPts fliplr(tPts)],...
            'YData',[(meanCvR{1}-stdCvR{1}) fliplr(meanCvR{1}+stdCvR{1})],...
            'FaceColor',[0 0.5 1],'FaceAlpha',0.4);
        uistack(p1,'bottom');
       
        xlim([tPts(1) tPts(end)]);
        ylim([-0.4 1]);
        line([tPts(1) tPts(end)],[0 0],'Color','k');
        line([0 0],[-0.4 1],'Color','k');
        xlabel('time (s)');
        ylabel('vR correlation');
        title('vR vs. DF/F correlation');
    end
end
 %% Save the plots
 
set(autoCPltPerTrial,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(autoCPltPerTrial,'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\autocorrPerTrial','-dpdf');

set(autoCPltAve,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(autoCPltAve, 'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\autocorrAve','-dpdf');

set(autoCPltPeak,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(autoCPltPeak, 'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\autocorrStats','-dpdf');

%% Plot the ROIs vs. the forward velocity

actScat = figure('units','normalized','outerposition',[0 0 1 1]);
actLines = figure('units','normalized','outerposition',[0 0 1 1]);
act3D = figure('units','normalized','outerposition',[0 0 1 1]);
actEx = figure('units','normalized','outerposition',[0 0 1 1]);

flyEx = 2;

vFEdges = linspace(0,1.5,16);
DFEdges = linspace(0,1.5,16);

for flyID=1:cond{1}.numFlies
    FlowDark = [];
    FlowCL = [];
    FlowFor = [];
    FlowRev = [];
    FlowRot = [];
    for trialID = 4:length(cond{2}.allFlyData{flyID}.Flow)
        % Find the translational open loop period
        OLPerTrans = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Trans > 0);
        
        % Find the forward translation
        ForPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Direction > 0);
        ForPer = intersect(OLPerTrans,ForPer);

        % Find the reverse translation
        RevPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Direction < 0);
        RevPer = intersect(OLPerTrans,RevPer);

        % Find the rotational open loop period
        OLPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.OLGain > 0);
        OLPerRot = setdiff(OLPer,OLPerTrans);

        % Find the clockwise rotations
        CWPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Direction > 0);
        CWPer = intersect(OLPerRot,CWPer);

        % Find the counter-clockwise rotations
        CCWPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Direction < 0);
        CCWPer = intersect(OLPerRot,CCWPer);
        
        % Find the dark periods
        DarkPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Direction == 0);
        
        % Find the closed loop period
        CLPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Closed > 0);
        CLPer = setdiff(CLPer,DarkPer);
               
        % Extract the imaging and velocity data
        ROIDat =  cond{2}.allFlyData{flyID}.Flow{trialID}.ROIaveMax-1;
        ROIDat(:,end) = [];
        ROIMax = max(ROIDat);
        vF = cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.vFFilt;
        OLGain = cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.OLGain;

        % Group data across trials
        movFor = find(vF > vFThresh);
        
        darkPlt = intersect(movFor,DarkPer);
        FlowDark = horzcat(FlowDark,vertcat(vF(darkPlt)',ROIMax(1,darkPlt)));
        
        CLPlt = intersect(movFor,CLPer);
        FlowCL = horzcat(FlowCL,vertcat(vF(CLPlt)',ROIMax(1,CLPlt)));
        
        ForPlt = intersect(movFor,ForPer);
        FlowFor = horzcat(FlowFor,vertcat(vF(ForPlt)',vertcat(OLGain(ForPlt)',ROIMax(1,ForPlt))));
        
        RevPlt = intersect(movFor,RevPer);
        FlowRev = horzcat(FlowRev,vertcat(vF(RevPlt)',vertcat(OLGain(RevPlt)',ROIMax(1,RevPlt))));
        
        RotPlt = intersect(movFor,OLPerRot);
        FlowRot = horzcat(FlowRot,vertcat(vF(RotPlt)',ROIMax(1,RotPlt)));
    end
    
    figure(actScat);
    
    subplot(5,cond{2}.numFlies,flyID);
    hold on;
    scatter(FlowDark(1,:),FlowDark(2,:),5,'k','filled');
    alpha(0.1);
    [darkFit darkFitError] = polyfit(FlowDark(1,:),FlowDark(2,:),1);
    linDark = polyval(darkFit,FlowDark(1,:),darkFitError);
    plot(FlowDark(1,:),linDark,'k');
    title('dark');
    xlim([0 1.5]);
    ylim([0 2]);
    text(0,2.5,strcat('fly #',num2str(flyID)));
    if flyID == 1
        ylabel('max DF/F');
    end
    
    subplot(5,cond{2}.numFlies,flyID+cond{2}.numFlies);
    hold on;
    scatter(FlowCL(1,:),FlowCL(2,:),5,'b','filled');
    alpha(0.1);
    [CLFit CLFitError] = polyfit(FlowCL(1,:),FlowCL(2,:),1);
    linCL = polyval(CLFit,FlowCL(1,:),CLFitError);
    plot(FlowCL(1,:),linCL,'b');
    title('CL');
    xlim([0 1.5]);
    ylim([0 2]);
    if flyID == 1
        ylabel('max DF/F');
    end
    
    subplot(5,cond{2}.numFlies,flyID+2*cond{2}.numFlies);
    hold on;
    scatter(FlowFor(1,:),FlowFor(3,:),5,'g','filled');
    alpha(0.1);
    [forFit forFitError] = polyfit(FlowFor(1,:),FlowFor(3,:),1);
    linFor = polyval(forFit,FlowFor(1,:),forFitError);
    plot(FlowFor(1,:),linFor,'g');
    title('for. flow');
    xlim([0 1.5]);
    ylim([0 2]);
    if flyID == 1
        ylabel('max DF/F');
    end
    
    subplot(5,cond{2}.numFlies,flyID+3*cond{2}.numFlies);
    hold on;
    scatter(FlowRev(1,:),FlowRev(3,:),5,'r','filled');
    alpha(0.1);
    [revFit revFitError] = polyfit(FlowRev(1,:),FlowRev(3,:),1);
    linRev = polyval(revFit,FlowRev(1,:),revFitError);
    plot(FlowRev(1,:),linRev,'r');
    title('rev. flow');
    xlim([0 1.5]);
    ylim([0 2]);
    if flyID == 1
        ylabel('max DF/F');
    end
    
    subplot(5,cond{2}.numFlies,flyID+4*cond{2}.numFlies);
    hold on;
    scatter(FlowRot(1,:),FlowRot(2,:),5,'m','filled');
    alpha(0.1);
    [rotFit rotFitError] = polyfit(FlowRot(1,:),FlowRot(2,:),1);
    linRot = polyval(rotFit,FlowRot(1,:),rotFitError);
    plot(FlowRot(1,:),linRot,'m');
    title('rot. flow');
    xlim([0 1.5]);
    ylim([0 2]);
    if flyID == 1
        ylabel('max DF/F');
    end
    xlabel('vF (cm/s)');
    
    figure(actLines);
    subplot(3,cond{2}.numFlies,flyID);
    hold on;
    h1 = histogram(FlowFor(1,:),vFEdges);
    h1.FaceColor = 'g';
    h1.FaceAlpha = 0.5;
    h1.EdgeColor = 'w';
    scatter(mean(FlowFor(1,:)),1000,'g','filled')
    line([mean(FlowFor(1,:))-std(FlowFor(1,:)) mean(FlowFor(1,:))+std(FlowFor(1,:))],...
        [1000 1000],'Color','g');
%     h2 = histogram(FlowCL(1,:),vFEdges);
%     h2.FaceColor = 'b';
%     h2.FaceAlpha = 0.2;
    scatter(mean(FlowCL(1,:)),980,'b','filled')
    line([mean(FlowCL(1,:))-std(FlowCL(1,:)) mean(FlowCL(1,:))+std(FlowCL(1,:))],...
        [980 980],'Color','b');
    h3 = histogram(FlowRev(1,:),vFEdges);
    h3.FaceColor = 'r';
    h3.FaceAlpha = 0.5;
    h3.EdgeColor = 'w';
    scatter(mean(FlowRev(1,:)),1000,'r','filled')
    line([mean(FlowRev(1,:))-std(FlowRev(1,:)) mean(FlowRev(1,:))+std(FlowRev(1,:))],...
        [1000 1000],'Color','r');
    xlim([0 1.5]);
    ylim([0 1000]);
    if flyID == 1
        ylabel('counts');
    end
    title(strcat('fly #',num2str(flyID)));
    
    subplot(3,cond{2}.numFlies,flyID+cond{2}.numFlies);
    hold on;
    ForColor = zeros(size(FlowFor));
    ForColor = ForColor + [0 0 0.5]';
    ForColor(2,:) = FlowFor(2,:)./max(FlowFor(2,:));
    RevColor = zeros(size(FlowRev));
    RevColor = RevColor + [0 0 0.5]';
    RevColor(1,:) = FlowRev(2,:)./max(FlowRev(2,:));
    scatter(FlowFor(1,:),FlowFor(3,:),2,'g','filled');
%     scatter(FlowCL(1,:),FlowCL(2,:),10,'b','filled');
    scatter(FlowRev(1,:),FlowRev(3,:),2,'r','filled');
    alpha(0.2);
    xlim([0 1.5]);
    ylim([0 2]);
    if flyID == 1
        ylabel('max DF/F');
    end
    
%     plot(FlowDark(1,:),linDark,'k');
%     plot(FlowCL(1,:),linCL,'b');
    plot(FlowFor(1,:),linFor,'g');
    plot(FlowRev(1,:),linRev,'r');
    xlabel('vF (cm/s)');
    
    subplot(3,cond{2}.numFlies,flyID+2*cond{2}.numFlies);
    hold on;
    h1 = histogram(FlowFor(3,:),DFEdges);
    h1.FaceColor = 'g';
    h1.FaceAlpha = 0.5;
    h1.EdgeColor = 'w';
    scatter(mean(FlowFor(3,:)),500,'g','filled')
    line([mean(FlowFor(3,:))-std(FlowFor(3,:)) mean(FlowFor(3,:))+std(FlowFor(3,:))],...
        [500 500],'Color','g');
%     h2 = histogram(FlowCL(3,:),DFEdges);
%     h2.FaceColor = 'b';
%     h2.FaceAlpha = 0.2;
    scatter(mean(FlowCL(2,:)),490,'b','filled')
    line([mean(FlowCL(2,:))-std(FlowCL(2,:)) mean(FlowCL(2,:))+std(FlowCL(2,:))],...
        [490 490],'Color','b');
    h3 = histogram(FlowRev(3,:),DFEdges);
    h3.FaceColor = 'r';
    h3.FaceAlpha = 0.5;
    h3.EdgeColor = 'w';
    scatter(mean(FlowRev(3,:)),500,'r','filled')
    line([mean(FlowRev(3,:))-std(FlowRev(3,:)) mean(FlowRev(3,:))+std(FlowRev(3,:))],...
        [500 500],'Color','r');
    xlim([0 1.5]);
    ylim([0 500]);
    if flyID == 1
        ylabel('counts');
    end
    xlabel('DF/F');
    
    figure(act3D)
        
    subplot(ceil(sqrt(cond{2}.numFlies)),ceil(sqrt(cond{2}.numFlies)), flyID);
    hold on;
    ForColor = zeros(size(FlowFor));
    ForColor = ForColor + [0 0 0.5]';
    ForColor(2,:) = FlowFor(2,:)./max(FlowFor(2,:));
    RevColor = zeros(size(FlowRev));
    RevColor = RevColor + [0 0 0.5]';
    RevColor(1,:) = FlowRev(2,:)./max(FlowRev(2,:));
    scatter3(FlowFor(1,:),FlowFor(2,:),FlowFor(3,:),2,ForColor');
    scatter3(FlowRev(1,:),-FlowRev(2,:),FlowRev(3,:),2,RevColor');
    alpha(0.2);
    xlim([0 1.5]);
    ylim([-100 100]);
    zlim([0 1.5]);
    view(-15,15);
    
    figure(actEx);
    if flyID == flyEx
        subplot(3,5,[11:12]);
        hold on;
        h1 = histogram(FlowFor(1,:),vFEdges);
        h1.FaceColor = 'g';
        h1.FaceAlpha = 0.5;
        h1.EdgeColor = 'w';
        scatter(mean(FlowFor(1,:)),1000,'g','filled');
        line([mean(FlowFor(1,:))-std(FlowFor(1,:)) mean(FlowFor(1,:))+std(FlowFor(1,:))],...
            [1000 1000],'Color','g');
        h3 = histogram(FlowRev(1,:),vFEdges);
        h3.FaceColor = 'm';
        h3.FaceAlpha = 0.5;
        h3.EdgeColor = 'w';
        scatter(mean(FlowRev(1,:)),1000,'m','filled');
        line([mean(FlowRev(1,:))-std(FlowRev(1,:)) mean(FlowRev(1,:))+std(FlowRev(1,:))],...
            [1000 1000],'Color','m');
        line([mean(FlowFor(1,:)) mean(FlowRev(1,:))],[1050 1050],'Color',[1 0.5 0]);
        xlim([0 1]);
        ylim([0 1050]);
        xlabel('vF (cm/s)');
        ylabel('counts');

        subplot(3,5,[1:2 6:7]);
        hold on;
%         plot(FlowFor(1,:),linFor,'g');
%         plot(FlowRev(1,:),linRev,'m');
        legend({'regressive flow','progressive flow'});
        legend('boxoff');
        scatter(FlowFor(1,:),FlowFor(3,:),20,'g','filled');
        scatter(FlowRev(1,:),FlowRev(3,:),20,'m','filled');
        alpha(0.2);
        xlim([0 1]);
        ylim([0 2]);
        ylabel('max DF/F');
        set(gca,'XTickLabel',[]);

        subplot(3,5,[3 8]);
        hold on;
        h1 = histogram(FlowFor(3,:),DFEdges,'Orientation','horizontal');
        h1.FaceColor = 'g';
        h1.FaceAlpha = 0.5;
        h1.EdgeColor = 'w';
        scatter(500,mean(FlowFor(3,:)),'g','filled')
        line([500 500],...
            [mean(FlowFor(3,:))-std(FlowFor(3,:)) mean(FlowFor(3,:))+std(FlowFor(3,:))],...
            [500 500],'Color','g');
        h3 = histogram(FlowRev(3,:),DFEdges,'Orientation','horizontal');
        h3.FaceColor = 'm';
        h3.FaceAlpha = 0.5;
        h3.EdgeColor = 'w';
        scatter(500,mean(FlowRev(3,:)),'m','filled')
        line([500 500],...
            [mean(FlowRev(3,:))-std(FlowRev(3,:)) mean(FlowRev(3,:))+std(FlowRev(3,:))],...
            'Color','m');
        line([525 525],[mean(FlowFor(3,:)) mean(FlowRev(3,:))],'Color','b');
        ylim([0 2]);
        xlim([0 525]);
        xlabel('counts');
        set(gca,'YTickLabel',[]);
    end
    
    subplot(3,5,[4:5 9:10 14:15]);
    hold on;
	scatter(1,mean(FlowFor(1,:))-mean(FlowRev(1,:)),100,[1.0 0.5 0],'filled');
    scatter(2,mean(FlowFor(3,:))-mean(FlowRev(3,:)),100,[0 0 1],'filled');
    line([1 2],[...
        mean(FlowFor(1,:))-mean(FlowRev(1,:)) ...
        mean(FlowFor(3,:))-mean(FlowRev(3,:))],...
        'Color',[0.8 0.8 0.8]);
    xlim([0 3]);
    ylim([-0.3 0.3]);
    line([0 3],[0 0],'Color','k');
    set(gca,'XTick',[1 2],'XTickLabel',{'mean vF diff.','mean DF/F diff.'});
    
end

% set(actScat,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(actScat,'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\ScatterConditions','-dpdf');
% 
% set(actLines,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(actLines,'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\ForVsRevFlowAll','-dpdf');

set(actEx,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(actEx,'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\ForVsRevExStats','-dpdf');

%% Compare the CL trajectory from the imaging to the trajectory calculated from the activity

innerPos = figure('units','normalized','outerposition',[0 0 1 1]);
for flyID=1:cond{2}.numFlies
    for trialID = 1:length(cond{2}.allFlyData{flyID}.Flow)
        
        % Find the dark periods
        DarkPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Direction == 0);
        
        % Find the closed loop period
        CLPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Closed > 0);
        CLPer = setdiff(CLPer,DarkPer);
        
        % Extract the imaging data
        FBDat =  cond{2}.allFlyData{flyID}.Flow{trialID}.ROIaveMax(:,CLPer)-1;
        maxF = max(FBDat);
        
        % Find the PVA
        num_ROIs = size(FBDat,1);
        angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
        angsraw = angsraw';
        clear meanAngRaw;
        clear meanIntRaw;
        for ts = 1:size(FBDat,2)
            meanAngRaw(ts) = circ_mean(angsraw,...
                squeeze(FBDat(:,ts)));
            meanIntRaw(ts) = circ_r(angsraw,...
                squeeze(FBDat(:,ts)));
        end
        
        % Find a fit between the DF max and the velocity
        tPts = cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.OffsetRotMatch(CLPer,1);
        vF = cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.vFFilt(CLPer);
        flyMov = find(vF>vFThresh);
        
        [p,S] = polyfit(vF(flyMov),maxF(flyMov)',1);

        % Plot the real and derived position
        XReal = cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.PosForMatch(CLPer);
        XReal = XReal - XReal(1);
        YReal = cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.PosLatMatch(CLPer);
        YReal = YReal - YReal(1);
        pReal = polyfit(XReal,YReal,1);
        XInner = zeros(length(vF),1);
        YInner = zeros(length(vF),1);
        tStep = mean(diff(tPts));
        for posStep = 2:length(XInner);
            XInner(posStep) = XInner(posStep-1)+tStep*cos(meanAngRaw(posStep-1))*(maxF(posStep-1)-p(2))./p(1);
            YInner(posStep) = YInner(posStep-1)+tStep*sin(meanAngRaw(posStep-1))*(maxF(posStep-1)-p(2))./p(1);
        end
        pInner = polyfit(XInner,YInner,1);
        u = XReal(end)+1i*YReal(end);
        v = XInner(end)+1i*YInner(end);
        angDiff = angle(u*conj(v));
        InnerRot = [cos(angDiff), -sin(angDiff);sin(angDiff), cos(angDiff)]*[XInner YInner]';
        
        subplot(6,cond{2}.numFlies,flyID+cond{2}.numFlies*(trialID-1));
        hold on;
        hold on;
        plot(XReal,YReal,'b');
        plot(InnerRot(1,:),InnerRot(2,:),'r');
        if flyID == 1 & trialID == 1
            legend({'real position','activity derived position'});
            legend('boxoff');
        end
        if flyID == 1
            ylabel(strcat('trial #',num2str(trialID)));
        else
            ylabel('y pos (cm)');
        end
        if trialID == 1
            title(strcat('fly #',num2str(flyID)));
        end
        xlabel('x pos (cm)');
        xlim([-20 20]);
        ylim([-20 20]);
        axis equal;
        
    end
end

set(innerPos,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(innerPos,'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\RealVsInnerPos2DVR','-dpdf');

%% Fit parameters for DF/F vs. vF relation


for flyID=1:cond{2}.numFlies
    
    innerPos = figure('units','normalized','outerposition',[0 0 1 1]);
    for trialID = 1:length(cond{2}.allFlyData{flyID}.Flow)
        
        % Find the dark periods
        DarkPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Direction == 0);
        
        % Find the closed loop period
        CLPer = find(cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.Closed > 0);
        CLPer = setdiff(CLPer,DarkPer);
        
        % Extract the imaging data
        FBDat =  cond{2}.allFlyData{flyID}.Flow{trialID}.ROIaveMax(:,CLPer)-1;
        maxF = max(FBDat);
        
        % Find the PVA
        num_ROIs = size(FBDat,1);
        angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
        angsraw = angsraw';
        clear meanAngRaw;
        clear meanIntRaw;
        for ts = 1:size(FBDat,2)
            meanAngRaw(ts) = circ_mean(angsraw,...
                squeeze(FBDat(:,ts)));
            meanIntRaw(ts) = circ_r(angsraw,...
                squeeze(FBDat(:,ts)));
        end
        
        % Find a fit between the DF max and the velocity
        tPts = cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.OffsetRotMatch(CLPer,1);
        vF = cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.vFFilt(CLPer);
        flyMov = find(vF>vFThresh);
        
        [p,S] = polyfit(vF(flyMov),maxF(flyMov)',1);

        % Plot the real and derived position
        XReal = cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.PosForMatch(CLPer);
        XReal = XReal - XReal(1);
        YReal = cond{2}.allFlyData{flyID}.Flow{trialID}.positionDatMatch.PosLatMatch(CLPer);
        YReal = YReal - YReal(1);
        pReal = polyfit(XReal,YReal,1);
        
        % Create an array of offsets and slopes to fit
        a = linspace(0,5,101);
        a(1) = [];
        b = linspace(0,1,21);
        b(1) = [];
        errorMat = zeros(length(a),length(b));
        tStep = mean(diff(tPts));
        
        for aTest = 1:length(a)
            for bTest = 1:length(b)
                
                XInner = zeros(length(vF),1);
                YInner = zeros(length(vF),1);
                
                for posStep = 2:length(XInner)
                    XInner(posStep) = XInner(posStep-1)+tStep*cos(meanAngRaw(posStep-1))*(maxF(posStep-1)-b(bTest))./a(aTest);
                    YInner(posStep) = YInner(posStep-1)+tStep*sin(meanAngRaw(posStep-1))*(maxF(posStep-1)-b(bTest))./a(aTest);
                end
                pInner = polyfit(XInner,YInner,1);
                u = XReal(end)+1i*YReal(end);
                v = XInner(end)+1i*YInner(end);
                angDiff = angle(u*conj(v));
                InnerRot = [cos(angDiff), -sin(angDiff);sin(angDiff), cos(angDiff)]*[XInner YInner]';
                
                errorMat(aTest,bTest) = sum(sqrt((InnerRot(1,:)'-XReal).^2+(InnerRot(2,:)'-YReal).^2));
            end
        end
        minInd = find(errorMat == min(min(errorMat)));
        bMin = ceil(minInd./length(a));
        aMin = minInd - (bMin-1)*length(a);
        
        XInner = zeros(length(vF),1);
        YInner = zeros(length(vF),1);

        for posStep = 2:length(XInner)
            XInner(posStep) = XInner(posStep-1)+tStep*cos(meanAngRaw(posStep-1))*(maxF(posStep-1)-b(bMin))./a(aMin);
            YInner(posStep) = YInner(posStep-1)+tStep*sin(meanAngRaw(posStep-1))*(maxF(posStep-1)-b(bMin))./a(aMin);
        end
        pInner = polyfit(XInner,YInner,1);
        u = XReal(end)+1i*YReal(end);
        v = XInner(end)+1i*YInner(end);
        angDiff = angle(u*conj(v));
        InnerRot = [cos(angDiff), -sin(angDiff);sin(angDiff), cos(angDiff)]*[XInner YInner]';
        
        subplot(3,3,1+mod((trialID-1),2)+3*floor((trialID-1)/2));
        hold on;
        hold on;
        plot(XReal,YReal,'b');
        plot(InnerRot(1,:),InnerRot(2,:),'r');
        if flyID == 1 & trialID == 1
            legend({'real position','activity derived position'});
            legend('boxoff');
        end
        if flyID == 1
            ylabel(strcat('trial #',num2str(trialID)));
        else
            ylabel('y pos (cm)');
        end
        if trialID == 1
            title(strcat('fly #',num2str(flyID)));
        end
        xlabel('x pos (cm)');
        axis equal;
        
        subplot(3,3,[3 6 9]);
        hold on;
        scatter(p(1),p(2),100,[0 0 trialID/6],'filled');
        scatter(a(aMin),b(bMin),100,[trialID/6 0 0] ,'filled');
        alpha(0.5);
        line([p(1) a(aMin)],[p(2) b(bMin)],'Color',[0.8 0.8 0.8]);
        xlim([0 max(a)]);
        xlabel('slope ((cm/s)/(DF/F)');
        ylabel('intersept (DF/F)');
        ylim([0 max(b)]);
        
    end
    
    set(innerPos,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(innerPos,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\InnerModelFit\InnerModelFit_Fly',num2str(flyID)),'-dpdf');
end


%% Plot 2D heading vs. forward velocity tuning curves for the individual FB ROIs

% Extract the relevant data
flyMap = {};
for flyID=1:cond{1}.numFlies
    actAll = [];
    vFAll = [];
    stripePosAll = [];
    for trialID = 1:length(cond{1}.allFlyData{flyID}.Stripe)
        
        % Find the dark periods
        DarkPer = find(cond{1}.allFlyData{flyID}.Stripe{trialID}.positionDatMatch.Direction == 0);
        
        % Find the closed loop period
        CLPer = find(cond{1}.allFlyData{flyID}.Stripe{trialID}.positionDatMatch.Closed > 0);
        CLPer = setdiff(CLPer,DarkPer);
        CLPer(end) = [];
        
        % Extract the appropriate data
        FBDat =  cond{1}.allFlyData{flyID}.Stripe{trialID}.ROIaveMax(:,CLPer)-1;
        vF = cond{1}.allFlyData{flyID}.Stripe{trialID}.positionDatMatch.vFFilt(CLPer);
        stripePos = cond{1}.allFlyData{flyID}.Stripe{trialID}.positionDatMatch.OffsetRotMatch(CLPer,2);
        
        actAll = horzcat(actAll,FBDat);
        vFAll = vertcat(vFAll,vF);
        stripePosAll = vertcat(stripePosAll,stripePos);
    end
    flyMap{flyID}.actAll = actAll;
    flyMap{flyID}.vFAll = vFAll;
    flyMap{flyID}.stripePosAll = stripePosAll;
end

% Create a 2D matrix and tuning curves out of the data
headingEdges = linspace(-pi,pi,17);
headingCents = headingEdges(1:end-1);
headingCents = headingCents+0.5*mean(diff(headingCents));
vFEdges = linspace(-0.25,2,10);
vFCents = vFEdges(1:end-1);
vFCents = vFCents+0.5*mean(diff(vFCents));

for flyID=1:cond{1}.numFlies
    % Initialize the objects and arrays
    mapMatrixAll = {};
    mapMatrixMean = zeros(8,length(headingEdges)-1,length(vFEdges)-1);
    headingTuneAll = {};
    headingTuneMean = zeros(8,length(headingEdges)-1);
    vFTuneAll = {};
    vFTuneMean = zeros(8,length(vFEdges)-1);
    
    % Bin the heading and velocity data
    headBin = discretize(flyMap{flyID}.stripePosAll,headingEdges);
    vFBin = discretize(flyMap{flyID}.vFAll,vFEdges);
    minN = 10;
    
    % Step through the ROIs
    for ROINow = 1:8
        % Populate the 2D array
        for headingStep = 1:length(headingEdges)-1
            for vFStep = 1:length(vFEdges)-1
                mapMatrixAll{ROINow}.heading{headingStep}.vF{vFStep} = [];
            end
        end
        for tPt = 1:length(flyMap{flyID}.vFAll)
            mapMatrixAll{ROINow}.heading{headBin(tPt)}.vF{vFBin(tPt)} = ...
                vertcat(mapMatrixAll{ROINow}.heading{headBin(tPt)}.vF{vFBin(tPt)},...
                flyMap{flyID}.actAll(ROINow,tPt));
        end
        for headingStep = 1:length(headingEdges)-1
            for vFStep = 1:length(vFEdges)-1
                mapMatrixMean(ROINow,headingStep,vFStep) = ...
                    mean(mapMatrixAll{ROINow}.heading{headingStep}.vF{vFStep});
            end
        end
        
        % Populate the linear arrays
        for headingStep = 1:length(headingEdges)-1
            headingTuneAll{ROINow}.heading{headingStep} = [];
        end
        for vFStep = 1:length(vFEdges)-1
            vFTuneAll{ROINow}.vF{vFStep} = [];
        end
        for tPt = 1:length(flyMap{flyID}.vFAll)
            headingTuneAll{ROINow}.heading{headBin(tPt)} = ...
                vertcat(headingTuneAll{ROINow}.heading{headBin(tPt)},...
                flyMap{flyID}.actAll(ROINow,tPt));
            vFTuneAll{ROINow}.vF{vFBin(tPt)} = ...
                vertcat(vFTuneAll{ROINow}.vF{vFBin(tPt)},...
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
        for vFStep = 1:length(vFEdges)-1
            if length(vFTuneAll{ROINow}.vF{vFStep}) > minN
                vFTuneMean(ROINow,vFStep) =...
                     mean(vFTuneAll{ROINow}.vF{vFStep});
            else
                vFTuneMean(ROINow,vFStep) = NaN;
            end
        end
    end
    flyMap{flyID}.mapMatrixAll = mapMatrixAll;
    flyMap{flyID}.mapMatrixMean = mapMatrixMean;
    flyMap{flyID}.headingTuneAll = headingTuneAll;
    flyMap{flyID}.headingTuneMean = headingTuneMean;
    flyMap{flyID}.vFTuneAll = vFTuneAll;
    flyMap{flyID}.vFTuneMean = vFTuneMean;
end

for flyID=1:cond{1}.numFlies
    vFHeadingMap = figure('units','normalized','outerposition',[0 0 1 1]);
    for ROINow = 1:8
        if ROINow <=4
            subplot(6,12,[1:2 13:14]+3*(ROINow-1));
            imagesc(vFCents,headingCents,flipud(squeeze(flyMap{flyID}.mapMatrixMean(ROINow,:,:))));
            axis off;
            title(strcat('ROI #',num2str(ROINow)));
            if ROINow == 1
                text(-0.25,-pi-0.5,strcat('fly #',num2str(flyID)),'FontSize',14);
            end
            
            subplot(6,12,[25:26]+3*(ROINow-1));
            plot(vFCents,flyMap{flyID}.vFTuneMean(ROINow,:),'Color',[1 0.5 0]);
            xlabel('vF (cm/s)');
            xlim([vFCents(1) vFCents(end)]);
            ylabel('DF/F');
            
            subplot(6,12,[3 15]+3*(ROINow-1));
            plot(flyMap{flyID}.headingTuneMean(ROINow,:),headingCents,'Color',[0 0 1]);
            xlabel('DF/F');            
            ylabel('heading (rad)');
            ylim([headingCents(1) headingCents(end)]);
            
        else
            subplot(6,12,[37:38 49:50]+3*(ROINow-5))
            imagesc(vFCents,headingCents,flipud(squeeze(flyMap{flyID}.mapMatrixMean(ROINow,:,:))));
            axis off;
            title(strcat('ROI #',num2str(ROINow)));
            
            subplot(6,12,[61:62]+3*(ROINow-5));
            plot(vFCents,flyMap{flyID}.vFTuneMean(ROINow,:),'Color',[1 0.5 0]);
            xlabel('vF (cm/s)');
            xlim([vFCents(1) vFCents(end)]);
            ylabel('DF/F');
            
            subplot(6,12,[39 51]+3*(ROINow-5));
            plot(flyMap{flyID}.headingTuneMean(ROINow,:),headingCents,'Color',[0 0 1]);
            xlabel('DF/F');            
            ylabel('heading (rad)');
            ylim([headingCents(1) headingCents(end)]);
        end
    end
    set(vFHeadingMap,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(vFHeadingMap,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\vFHeadingMap\vFHeadingMap_fly',num2str(flyID)),'-dpdf');
    delete(vFHeadingMap);

end

%% Find off and on time constants that lead to the highest correlation values

% Specify processing parameters
GtOnAll = [0.001:0.001:0.1];
GtOffAll = [0.01:0.01:0.5];

condID = 1;

sgolayOrder = 3;
sgolayFrames = 11;

allGtOn = zeros(cond{1}.numFlies,3,3);
allGtOff = zeros(cond{1}.numFlies,3,3);
bestCC = zeros(cond{1}.numFlies,3,3);

for flyID = 1:cond{condID}.numFlies
    for trialType = 1:3
        if trialType == 1
            trialName = 'StripeWGround';
        elseif trialType == 2
            trialName = 'MSClutterWGround';
        else
            trialName = 'MSClutterWGroundWCyl';
        end
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

            % Find the PVA, the forward velocity, and the max DF
            tAll = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll - tAll(1);
            tStep = mean(diff(tAll));

            FBDat =  cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;
            maxDF = max(FBDat);
            maxDF = sgolayfilt(maxDF,sgolayOrder,sgolayFrames);

            vF = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vF;
            
            allCCs = zeros(length(GtOnAll),length(GtOffAll));
            for onTest = 1:length(GtOnAll)
                GtOn = GtOnAll(onTest)/log(2);
                for offTest = 1:length(GtOffAll)
                    GtOff = GtOffAll(offTest)/log(2);
                    vFConv = conv(vF',abs(exp(-tAll/GtOn)-exp(-tAll/GtOff)));
                    CCs = corrcoef(maxDF(2:end),vFConv(1:length(vF)));
                    allCCs(onTest,offTest) = CCs(2,1);
                end
            end

            [maxCC maxInd]=max(allCCs(:));
            [X Y]=ind2sub(size(allCCs),maxInd);
            allGtOn(flyID,trialType,trialID) = GtOnAll(X);
            allGtOff(flyID,trialType,trialID) = GtOffAll(Y);
            bestCC(flyID,trialType,trialID) = maxCC;
        end
    end
end

bestGtOn = median(allGtOn(:));
bestGtOff = median(allGtOff(:));

%% Plot 2D heading vs. forward velocity tuning curves for the individual FB ROIs - for the 2D stims and conv. vFs

minN = 10;

% Extract the relevant data
flyMap = {};
for flyID=1:cond{1}.numFlies
    for trialType = 1:3
        if trialType == 1
            trialName = 'StripeWGround';
        elseif trialType == 2
            trialName = 'MSClutterWGround';
        else
            trialName = 'MSClutterWGroundWCyl';
        end
        for trialID = 1:length(cond{1}.allFlyData{flyID}.(trialName))

            % Extract the appropriate data
            FBDat =  cond{1}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;
            tAll = cond{1}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll - tAll(1);
            vF = cond{1}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vF;
            vF = conv(vF',abs(exp(-tAll/bestGtOn)-exp(-tAll/bestGtOff)));
            vF(size(FBDat,2):end) = [];
            heading = cond{1}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,2);
            
            flyMap{flyID}.(trialName){trialID}.FBDat = FBDat;
            flyMap{flyID}.(trialName){trialID}.vF = vF;
            flyMap{flyID}.(trialName){trialID}.heading = heading;
        end
    end
end

% Create a 2D matrix and tuning curves out of the data
headingEdges = linspace(-pi,pi,17);
headingCents = headingEdges(1:end-1);
headingCents = headingCents+0.5*mean(diff(headingCents));
vFEdges = linspace(-0.25,6,26);
vFCents = vFEdges(1:end-1);
vFCents = vFCents+0.5*mean(diff(vFCents));

for flyID=1:cond{1}.numFlies
    for trialType = 1:3
        if trialType == 1
            trialName = 'StripeWGround';
        elseif trialType == 2
            trialName = 'MSClutterWGround';
        else
            trialName = 'MSClutterWGroundWCyl';
        end
        for trialID = 1:length(cond{1}.allFlyData{flyID}.(trialName))
            % Initialize the objects and arrays
            mapMatrixAll = {};
            mapMatrixMean = [];
            vFTuneAll = {};
            vFTuneMean = [];
            
            % Bin the heading and velocity data
            headBin = discretize(flyMap{flyID}.(trialName){trialID}.heading,headingEdges);
            vFBin = discretize(flyMap{flyID}.(trialName){trialID}.vF,vFEdges);

            % Step through the ROIs
            for ROINow = 1:8
                % Populate the 2D array
                for headingStep = 1:length(headingEdges)-1
                    for vFStep = 1:length(vFEdges)-1
                        mapMatrixAll{ROINow}.heading{headingStep}.vF{vFStep} = [];
                    end
                end
                for tPt = 1:length(flyMap{flyID}.(trialName){trialID}.vF)
                    mapMatrixAll{ROINow}.heading{headBin(tPt)}.vF{vFBin(tPt)} = ...
                        vertcat(mapMatrixAll{ROINow}.heading{headBin(tPt)}.vF{vFBin(tPt)},...
                        flyMap{flyID}.(trialName){trialID}.FBDat(ROINow,tPt));
                end
                for headingStep = 1:length(headingEdges)-1
                    for vFStep = 1:length(vFEdges)-1
                        mapMatrixMean(ROINow,headingStep,vFStep) = ...
                            mean(mapMatrixAll{ROINow}.heading{headingStep}.vF{vFStep});
                    end
                end
                
                % Populate the linear arrays
                for vFStep = 1:length(vFEdges)-1
                    vFTuneAll{ROINow}.vF{vFStep} = [];
                end
                for tPt = 1:length(flyMap{flyID}.(trialName){trialID}.vF)
                    vFTuneAll{ROINow}.vF{vFBin(tPt)} = ...
                        vertcat(vFTuneAll{ROINow}.vF{vFBin(tPt)},...
                        flyMap{flyID}.(trialName){trialID}.FBDat(ROINow,tPt));
                end
                for vFStep = 1:length(vFEdges)-1
                    if length(vFTuneAll{ROINow}.vF{vFStep}) > minN
                        vFTuneMean(ROINow,vFStep) =...
                             mean(vFTuneAll{ROINow}.vF{vFStep});
                    else
                        vFTuneMean(ROINow,vFStep) = NaN;
                    end
                end
            end
            flyMap{flyID}.(trialName){trialID}.mapMatrixAll = mapMatrixAll;
            flyMap{flyID}.(trialName){trialID}.mapMatrixMean = mapMatrixMean;
            flyMap{flyID}.(trialName){trialID}.vFTuneAll = vFTuneAll;
            flyMap{flyID}.(trialName){trialID}.vFTuneMean = vFTuneMean;
        end
    end
end

trialName = 'StripeWGround';

for flyID=1:cond{1}.numFlies    
    vFHeadingMap = figure('units','normalized','outerposition',[0 0 1 1]);
    for trialID = 1:length(cond{1}.allFlyData{flyID}.(trialName))
        for ROINow = 1:8

            subplot(3,8,ROINow+8*(trialID-1));
            imagesc(vFCents,headingCents,flipud(squeeze(flyMap{flyID}.(trialName){trialID}.mapMatrixMean(ROINow,:,:))));
            title(strcat('ROI #',num2str(ROINow)));
            if ROINow == 1
                tLab = text(-1.5,-1.1*pi,strcat('trial #',num2str(trialID)));
                if trialID == 1
                    text(-2,-1.5*pi,strcat('fly #',num2str(flyID)),'FontSize',14);
                end
            end
            xlim([0 2]);
            if trialID == 3
                xlabel('vF (cm/s)');
            end
            if ROINow == 1
                ylabel('heading (rad)');
            end

        end
    end
    set(vFHeadingMap,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(vFHeadingMap,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\vFHeadingMap2DStripe_fly',num2str(flyID)),'-dpdf');
end



% vFTuning = figure('units','normalized','outerposition',[0 0 1 1]);
% 
% for flyID=1:cond{1}.numFlies    
%     for trialID = 1:length(cond{1}.allFlyData{flyID}.(trialName))
%         for ROINow = 1:8
% 
%             subplot(6,8,ROINow+8*(trialID-1)+24*(flyID-1));
%             plot(vFCents,squeeze(flyMap{flyID}.(trialName){trialID}.vFTuneMean(ROINow,:)));
%            
%         end
%     end
% end

%% Plot the trajectories color coded by the offsets

condID = 1;

offsetEdges = linspace(-pi,pi,33);
offsetCents = offsetEdges(1:end-1);
offsetCents = offsetCents+0.5*mean(diff(offsetCents));
offsetColors = zeros(length(offsetCents),3);
colorSpan = zeros(length(offsetCents),1);
colorSpan(1:floor(length(offsetCents)/2)) = linspace(0,1,floor(length(offsetCents)/2));
colorSpan(floor(length(offsetCents)/2)+1:end) = linspace(1,0,length(offsetCents)-floor(length(offsetCents)/2));
offsetColors(:,1) = colorSpan;
offsetColors(:,2) = circshift(colorSpan,floor(length(offsetCents)/3));
offsetColors(:,3) = circshift(colorSpan,floor(2*length(offsetCents)/3));

for flyID=1:cond{condID}.numFlies
    offsetCoding = figure('units','normalized','outerposition',[0 0 1 1]);
    for trialType = 1:3
        if trialType == 1
            trialName = 'StripeWGround';
        elseif trialType == 2
            trialName = 'MSClutterWGround';
        else
            trialName = 'MSClutterWGroundWCyl';
        end
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

            % Extract the appropriate data
            FBDat =  cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;
            forPos = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetForMatch;
            latPos = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetLatMatch;
            vF = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vF;
            vF = conv(vF',abs(exp(-tAll/bestGtOn)-exp(-tAll/bestGtOff)));
            vF(size(FBDat,2):end) = [];
            heading = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,2);
            
            % Find the PVA
            num_ROIs = size(FBDat,1);
            angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
            angsraw = angsraw';
            clear PVA;
            clear PVAStren;
            for ts = 1:size(FBDat,2)
                PVA(ts) = circ_mean(angsraw,...
                    squeeze(FBDat(:,ts)));
                PVAStren(ts) = circ_r(angsraw,...
                    squeeze(FBDat(:,ts)));
            end
            
            offset = mod(PVA'-heading,2*pi)-pi;
            offsetIDs = discretize(offset,offsetEdges);
            pltColors = zeros(length(offset),3);
            for pltPt = 1:length(offset)
                pltColors(pltPt,:) = offsetColors(offsetIDs(pltPt),:);
            end
            
            
            % Plot the trajectories colored by the offset
            subplot(3,4,1+trialType+4*(trialID-1));
            scatter(forPos,latPos,4,pltColors,'filled');
            rectangle('Position',[-75 -75 150 150],'Curvature',1,'EdgeColor','k');
            if trialType == 1
                rectangle('Position',[75 0 2 10],'FaceColor','c');
            elseif trialType == 2
                rectangle('Position',[-75 -75 150 150],'Curvature',1,'EdgeColor','c','LineStyle','--');
            elseif trialType == 3
                rectangle('Position',[-0.5 -0.5 1 1],'Curvature',1,'FaceColor','c');
            end
                
            axis equal;
            if flyID == 1
                xlim([-25 25]);
            end
            axis off;
            
            
            if trialID == 1
                title(trialName);
            end
            if trialType == 1
                xlimNow = get(gca,'XLim');
                ylimNow = get(gca,'YLim');
                text(xlimNow(1),ylimNow(2),strcat('trial #',num2str(trialID)));
                if trialID == 1
                    line([xlimNow(1) xlimNow(1)+10],[ylimNow(1) ylimNow(1)],'Color','k');
                    text(xlimNow(1),ylimNow(1)+5,'10 cm');
                end
            end
            
        end
    end
    subplot(3,4,5);
    rectangle('Position',[-1 -1 2 2],'Curvature',1,'EdgeColor','k');
    scatter(cos(offsetCents),sin(offsetCents),100,offsetColors,'filled');
    title('offset');
    axis equal;
    axis off;
    
    set(offsetCoding,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(offsetCoding,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\OffsetCodedTraj_Fly#',num2str(flyID)),'-dpdf');
end

%% Plot the trajectories color coded by the offsets - stripe & clutter

condID = 1;

offsetEdges = linspace(-pi,pi,33);
offsetCents = offsetEdges(1:end-1);
offsetCents = offsetCents+0.5*mean(diff(offsetCents));
offsetColors = zeros(length(offsetCents),3);
colorSpan = zeros(length(offsetCents),1);
colorSpan(1:floor(length(offsetCents)/2)) = linspace(0,1,floor(length(offsetCents)/2));
colorSpan(floor(length(offsetCents)/2)+1:end) = linspace(1,0,length(offsetCents)-floor(length(offsetCents)/2));
offsetColors(:,1) = colorSpan;
offsetColors(:,2) = circshift(colorSpan,floor(length(offsetCents)/3));
offsetColors(:,3) = circshift(colorSpan,floor(2*length(offsetCents)/3));

offsetCodingAll = figure('units','normalized','outerposition',[0 0 1 1]);
offsetCodingClutter = figure('units','normalized','outerposition',[0 0 1 1]);

for flyID=1:cond{condID}.numFlies    

    for trialType = 1:2
        if trialType == 1
            trialName = 'StripeWGround';
        elseif trialType == 2
            trialName = 'MSClutterWGround';
        end
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

            % Extract the appropriate data
            FBDat =  cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;
            forPos = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetForMatch;
            latPos = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetLatMatch;
            vF = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vF;
            vF = conv(vF',abs(exp(-tAll/bestGtOn)-exp(-tAll/bestGtOff)));
            vF(size(FBDat,2):end) = [];
            heading = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,2);

            % Find the PVA
            num_ROIs = size(FBDat,1);
            angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
            angsraw = angsraw';
            clear PVA;
            clear PVAStren;
            for ts = 1:size(FBDat,2)
                PVA(ts) = circ_mean(angsraw,...
                    squeeze(FBDat(:,ts)));
                PVAStren(ts) = circ_r(angsraw,...
                    squeeze(FBDat(:,ts)));
            end

            offset = mod(PVA'-heading,2*pi)-pi;
            offsetIDs = discretize(offset,offsetEdges);
            pltColors = zeros(length(offset),3);
            for pltPt = 1:length(offset)
                pltColors(pltPt,:) = offsetColors(offsetIDs(pltPt),:);
            end


            figure(offsetCodingAll);
            % Plot the trajectories colored by the offset
            subplot(4,6,[1+flyID+12*(trialType-1) 7+flyID+12*(trialType-1)]);
            hold on;
            scatter(forPos,latPos,3,pltColors,'filled');
            rectangle('Position',[-75 -75 150 150],'Curvature',1,'EdgeColor','k');
            if trialType == 1
                rectangle('Position',[75 0 2 10],'FaceColor','c');
                line([-5 -70*cos(pi/3)],[0 -70*sin(pi/3)],'Color','k');
                line([-5 -70*cos(pi/3)],[0 70*sin(pi/3)],'Color','k');
            elseif trialType == 2
                rectangle('Position',[-75 -75 150 150],'Curvature',1,'EdgeColor','c','LineStyle','--');
            end
            axis equal;
            axis off;


            if flyID == 1
                title(trialName);
            end
            if trialType == 1 & trialID == 1 & flyID == 1
                line([30 40],[-40 -40],'Color','k');
                text(30,-35,'10 cm');
            end
            
            if trialType == 2
                figure(offsetCodingClutter);
                % Plot the trajectories colored by the offset
                subplot(floor(sqrt(cond{condID}.numFlies)),ceil(sqrt(cond{condID}.numFlies)),flyID);
                hold on;
                scatter(forPos,latPos,3,pltColors,'filled');
                scatter(-5,0,100,'k','p');
                axis equal;
                axis off;
            end
        end
    end
end
figure(offsetCodingAll);
subplot(4,6,[7 13]);
rectangle('Position',[-1 -1 2 2],'Curvature',1,'EdgeColor','k');
scatter(cos(offsetCents),sin(offsetCents),100,offsetColors,'filled');
title('offset');
axis equal;
axis off;
    
set(offsetCodingAll,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(offsetCodingAll,'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\OffsetCodedTrajAll','-dpdf');

set(offsetCodingClutter,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(offsetCodingClutter,'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\OffsetCodedTrajClutter','-dpdf');

%% Plot the trajectories color coded by the offsets - clutter w/cylinder

condID = 1;

offsetEdges = linspace(-pi,pi,33);
offsetCents = offsetEdges(1:end-1);
offsetCents = offsetCents+0.5*mean(diff(offsetCents));
offsetColors = zeros(length(offsetCents),3);
colorSpan = zeros(length(offsetCents),1);
colorSpan(1:floor(length(offsetCents)/2)) = linspace(0,1,floor(length(offsetCents)/2));
colorSpan(floor(length(offsetCents)/2)+1:end) = linspace(1,0,length(offsetCents)-floor(length(offsetCents)/2));
offsetColors(:,1) = colorSpan;
offsetColors(:,2) = circshift(colorSpan,floor(length(offsetCents)/3));
offsetColors(:,3) = circshift(colorSpan,floor(2*length(offsetCents)/3));

offsetCodingCylinder = figure('units','normalized','outerposition',[0 0 1 1]);
trialName = 'MSClutterWGroundWCyl';

for flyID=1:cond{condID}.numFlies    

    
    for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

        % Extract the appropriate data
        FBDat =  cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;
        forPos = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetForMatch;
        latPos = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetLatMatch;
        vF = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vF;
        vF = conv(vF',abs(exp(-tAll/bestGtOn)-exp(-tAll/bestGtOff)));
        vF(size(FBDat,2):end) = [];
        heading = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,2);
        
        % Calculate the object heading
        cylHeading = atan2(latPos,-forPos);
        totHeading = mod(heading+cylHeading+pi,2*pi)-pi;

        % Find the PVA
        num_ROIs = size(FBDat,1);
        angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
        angsraw = angsraw';
        clear PVA;
        clear PVAStren;
        for ts = 1:size(FBDat,2)
            PVA(ts) = circ_mean(angsraw,...
                squeeze(FBDat(:,ts)));
            PVAStren(ts) = circ_r(angsraw,...
                squeeze(FBDat(:,ts)));
        end

        offset = mod(PVA'-heading,2*pi)-pi;
        offsetCyl = mod(PVA'-totHeading,2*pi)-pi;
        offsetIDs = discretize(offset,offsetEdges);
        offsetIDsCyl = discretize(offsetCyl,offsetEdges);
        pltColors = zeros(length(offset),3);
        pltColorsCyl = zeros(length(offsetCyl),3);
        for pltPt = 1:length(offset)
            pltColors(pltPt,:) = offsetColors(offsetIDs(pltPt),:);
            pltColorsCyl(pltPt,:) = offsetColors(offsetIDsCyl(pltPt),:);
        end

        % Plot the trajectories colored by the offset
        subplot(4,1+cond{condID}.numFlies,[1+flyID 2+cond{condID}.numFlies+flyID]);
        hold on;
        scatter(forPos,latPos,3,pltColors,'filled');
        rectangle('Position',[-75 -75 150 150],'Curvature',1,'EdgeColor','k');
        rectangle('Position',[-0.5 -0.5 1 1],'Curvature',1,'FaceColor','c');
        axis equal;
        axis off;
        xlim([-10 10]);
        ylim([-10 10]);
        if trialID == 1 & flyID == 1
            line([5 10],[-10 -10],'Color','k');
            text(7,-8,'10 cm');
            text(-15,10,'global heading');
        end
        
        subplot(4,1+cond{condID}.numFlies,[3+2*cond{condID}.numFlies+flyID 4+3*cond{condID}.numFlies+flyID]);
        hold on;
        scatter(forPos,latPos,3,pltColorsCyl,'filled');
        rectangle('Position',[-75 -75 150 150],'Curvature',1,'EdgeColor','k');
        rectangle('Position',[-0.5 -0.5 1 1],'Curvature',1,'FaceColor','c');
        axis equal;
        axis off;
        xlim([-10 10]);
        ylim([-10 10]);
        if trialID == 1 & flyID == 1
            text(-15,10,'local heading');
        end

    end
end

subplot(4,1+cond{condID}.numFlies,[2+cond{condID}.numFlies 3+2*cond{condID}.numFlies]);
rectangle('Position',[-1 -1 2 2],'Curvature',1,'EdgeColor','k');
scatter(cos(offsetCents),sin(offsetCents),100,offsetColors,'filled');
title('offset');
axis equal;
axis off;

set(offsetCodingCylinder,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(offsetCodingCylinder,'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\OffsetCodedTrajCyl','-dpdf');

%% Bin the offsets by heading - clutter

condID = 1;

PVAThresh = 0.2;

offsetEdges = linspace(-pi,pi,9);
offsetCents = offsetEdges(1:end-1);
offsetCents = offsetCents+0.5*mean(diff(offsetCents));
offsetColors = zeros(length(offsetCents),3);
colorSpan = zeros(length(offsetCents),1);
colorSpan(1:floor(length(offsetCents)/2)) = linspace(0,1,floor(length(offsetCents)/2));
colorSpan(floor(length(offsetCents)/2)+1:end) = linspace(1,0,length(offsetCents)-floor(length(offsetCents)/2));
offsetColors(:,1) = colorSpan;
offsetColors(:,2) = circshift(colorSpan,floor(length(offsetCents)/3));
offsetColors(:,3) = circshift(colorSpan,floor(2*length(offsetCents)/3));

for flyID=1:cond{condID}.numFlies    

    trialName = 'MSClutterWGround';
    offsetCodingHist = figure('units','normalized','outerposition',[0 0 1 1]);
    
    NAll = zeros(length(offsetCents),length(offsetCents));
    
    for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

        % Extract the appropriate data
        FBDat =  cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;
        forPos = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetForMatch;
        latPos = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetLatMatch;
        heading = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,2);

        % Find the PVA
        num_ROIs = size(FBDat,1);
        angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
        angsraw = angsraw';
        clear PVA;
        clear PVAStren;
        for ts = 1:size(FBDat,2)
            PVA(ts) = circ_mean(angsraw,...
                squeeze(FBDat(:,ts)));
            PVAStren(ts) = circ_r(angsraw,...
                squeeze(FBDat(:,ts)));
        end

        headingIDs = discretize(heading,offsetEdges);
        
        offset = mod(PVA'-heading,2*pi)-pi;
        offsetIDs = discretize(offset,offsetEdges);
        pltColors = zeros(length(offset),3);
        for pltPt = 1:length(offset)
            pltColors(pltPt,:) = offsetColors(offsetIDs(pltPt),:);
        end

        % Plot the trajectories colored by the offset
        for angNow = 1:length(offsetCents)
            PVANow = intersect(find(headingIDs == angNow),find(PVAStren > PVAThresh));
            
            if ~isnan(PVANow)
                subplot(4,9,1+angNow+9*(trialID-1));
                hold on;
                quiver(0,0,cos(offsetCents(angNow)),sin(offsetCents(angNow)),'Color','k');
                [NPVA,edgesPVA] = histcounts(offset(PVANow),offsetEdges);
                NAll(angNow,:) = NAll(angNow,:)+NPVA;
                text(-1.5,1.5,strcat('Nmax = ',num2str(max(NPVA))));
                NPVA = NPVA./max(NPVA);
                for histBin = 1:length(edgesPVA)-1
                    angPlt = linspace(edgesPVA(histBin),edgesPVA(histBin+1),5);
                    radPtsX = cos(angPlt);
                    radPtsY = sin(angPlt);
                    patch([radPtsX 1.2*fliplr(radPtsX)],[radPtsY 1.2*fliplr(radPtsY)],[1-NPVA(histBin) 1-NPVA(histBin) 1-NPVA(histBin)],...
                        'EdgeColor',offsetColors(histBin,:));
                end
                axis equal;
                axis off;
            end
        end
        
        subplot(4,9,1+9*(trialID-1));
        hold on;
        scatter(forPos,latPos,3,pltColors,'filled');
        scatter(-5,0,25,'k','p');
        axis equal;
        axis off;
    end
    
    for angNow = 1:length(offsetCents)
        subplot(4,9,28+angNow);
        hold on;
        quiver(0,0,cos(offsetCents(angNow)),sin(offsetCents(angNow)),'Color','k');
        NNow = NAll(angNow,:);
        NNow = NNow./max(NNow);
        
        for histBin = 1:length(edgesPVA)-1
            angPlt = linspace(edgesPVA(histBin),edgesPVA(histBin+1),5);
            radPtsX = cos(angPlt);
            radPtsY = sin(angPlt);
            patch([radPtsX 1.2*fliplr(radPtsX)],[radPtsY 1.2*fliplr(radPtsY)],[1-NNow(histBin) 1-NNow(histBin) 1-NNow(histBin)],...
                'EdgeColor',offsetColors(histBin,:));
        end
        axis equal;
        axis off;
    end
end

    
% set(offsetCodingAll,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(offsetCodingAll,'C:\Users\turnerevansd\Documents\RawAnalysis\37G12\Figures\OffsetCodedTrajAll','-dpdf');
