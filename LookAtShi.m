%% Clear out old data
clear;
clc;

cond = FlyDatLoad(1);

%% Plot an example for each condition

% Find the ROI angles
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));

% Generate figures
compFigDark = figure('units','normalized','outerposition',[0 0 1 1]);
compFigCL = figure('units','normalized','outerposition',[0 0 1 1]);
compFigOL = figure('units','normalized','outerposition',[0 0 1 1]);

% Step through the conditions
for condID = 1:length(cond)
    
    if condID == 1
        flyID = 3;
        trialID = 2;
    elseif condID == 2
        flyID = 4;
        trialID = 2;
    else
        flyID = 1;
        trialID = 3;
    end
    
    % Extract the behavioral data
    tAll = cond{condID}.allFlyData{flyID}.All_30C{trialID}.positionDatMatch.OffsetRotMatch(:,1);
    tAll = tAll-tAll(1);
    heading = -pi/180*cond{condID}.allFlyData{flyID}.All_30C{trialID}.positionDatMatch.PosRotMatch;
    stripe = -cond{condID}.allFlyData{flyID}.All_30C{trialID}.positionDatMatch.OffsetRotMatch(:,2);
        
    % Find the dark periods
    DarkPer = find(cond{condID}.allFlyData{flyID}.All_30C{trialID}.positionDatMatch.Direction == 0);

    % Find the closed loop period
    CLPer = find(cond{condID}.allFlyData{flyID}.All_30C{trialID}.positionDatMatch.Closed > 0);
    CLPer = setdiff(CLPer,DarkPer);
    CLPer(end) = [];

    % Find the rotational open loop period
    OLPer = find(cond{condID}.allFlyData{flyID}.All_30C{trialID}.positionDatMatch.OLGain > 0);

    % Step through the visual conditions
    for stim = 1:3
        if stim == 1
            perNow = DarkPer;
            figure(compFigDark);
        elseif stim == 2
            perNow = CLPer;
            figure(compFigCL);
        else
            perNow = OLPer;
            figure(compFigOL);
        end
            
        subplot(3,1,condID);
        hold on;
        imagesc(tAll(perNow),angs,cond{condID}.allFlyData{flyID}.All_30C{trialID}.ROIaveMax(:,perNow)-1);
        plot(tAll(perNow),-heading(perNow),'k');
        if stim > 1
            plot(tAll(perNow),-stripe(perNow),'b');
        end
        xlim([tAll(perNow(1)) tAll(perNow(end))]);
        ylim([-pi pi]);
        caxis([0 1]);
        colorbar;
        title(cond{condID}.name);
        if condID == 1
            legend({'heading','stripe position'});
            legend('boxoff');
        end
        if condID == 3
            xlabel('time');
        end
        ylabel('position');
        
        colormap(brewermap(64, 'Blues'));
    end
end

set(compFigDark,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(compFigDark,'C:\Users\turnerevansd\Documents\RawAnalysis\Shi\Figures\ShiExample_Dark','-dpdf');
set(compFigCL,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(compFigCL,'C:\Users\turnerevansd\Documents\RawAnalysis\Shi\Figures\ShiExample_CL','-dpdf');
% set(compFigOL,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(compFigOL,'C:\Users\turnerevansd\Documents\RawAnalysis\Shi\Figures\ShiExample_OL','-dpdf');

%% Sort the activity by the visual conditions

allAct = {};

heatCond = 'All_30C';

for condID = 1:length(cond)
   allAct{condID}.name = cond{condID}.name;
   
   % Step through the flies
   for flyID = 1:cond{condID}.numFlies
       
      % Sort the data across trials
      for trialID = 1:length(cond{condID}.allFlyData{flyID}.(heatCond))
         
         % Find the dark period
         darkPer = find(cond{condID}.allFlyData{flyID}.(heatCond){trialID}.positionDatMatch.Direction == 0);
         
         % Find the open loop period
         OLPer = find(cond{condID}.allFlyData{flyID}.(heatCond){trialID}.positionDatMatch.OLGain > 0);
         
         % Find the clockwise rotations
         CWPer = find(cond{condID}.allFlyData{flyID}.(heatCond){trialID}.positionDatMatch.Direction > 0);
         CWPer = intersect(OLPer,CWPer);
         % Find the counter-clockwise rotations
         CCWPer = find(cond{condID}.allFlyData{flyID}.(heatCond){trialID}.positionDatMatch.Direction < 0);
         CCWPer = intersect(OLPer,CCWPer);
         
         % Find the closed loop period
         CLPer = find(cond{condID}.allFlyData{flyID}.(heatCond){trialID}.positionDatMatch.Direction > 0);
         CLPer = setdiff(CLPer,CWPer);
         
         % Extract the behavioral parameters
         vR = cond{condID}.allFlyData{flyID}.(heatCond){trialID}.positionDatMatch.vRot;
         vF = cond{condID}.allFlyData{flyID}.(heatCond){trialID}.positionDatMatch.vF;
         tPts = cond{condID}.allFlyData{flyID}.(heatCond){trialID}.positionDatMatch.OffsetRotMatch(:,1);
         stripePos = cond{condID}.allFlyData{flyID}.(heatCond){trialID}.positionDatMatch.OffsetRotMatch(:,2);
         heading = cond{condID}.allFlyData{flyID}.(heatCond){trialID}.positionDatMatch.PosRotMatch;
         
         % Get rid of the last time point so that I have velocities
         % throughout
         if CWPer(end) > length(vR)
             CWPer(end) = [];
         elseif CCWPer(end) > length(vR)
             CCWPer(end) = [];
         elseif darkPer(end) > length(vR)
             darkPer(find(darkPer > OLPer(end))) = [];
         end
         
         % Sort the data by period and color
         for periodID = 1:4
             
                 if periodID == 1
                     perNow = darkPer;
                     allAct{condID}.fly{flyID}.period{periodID}.type = 'dark';
                     allAct{condID}.fly{flyID}.period{periodID}.color = [0 0 0];
                 elseif periodID == 2
                     perNow = CLPer;
                     allAct{condID}.fly{flyID}.period{periodID}.type = 'CL';
                     allAct{condID}.fly{flyID}.period{periodID}.color = [0 0 1];
                 elseif periodID == 3
                     perNow = CWPer;
                     allAct{condID}.fly{flyID}.period{periodID}.type = 'CW';
                     allAct{condID}.fly{flyID}.period{periodID}.color = [0.5 0 1];
                 elseif periodID == 4
                     perNow = CCWPer;
                     allAct{condID}.fly{flyID}.period{periodID}.type = 'CCW';
                     allAct{condID}.fly{flyID}.period{periodID}.color = [0 0.5 1];
                 end
                 allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vR = vR(perNow);
                 allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vF = vF(perNow);
                 
                 % Find the time points, stripe position, and heading
                 allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.tPts = tPts(perNow);
                 allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.stripePos = stripePos(perNow);
                 allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.heading = heading(perNow);
             % Find the activity
             allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.act = ...
                 cond{condID}.allFlyData{flyID}.(heatCond){trialID}.ROIaveMax(:,perNow);
         end
      end
   end
end

%% Plot the PVA vs. position for the dark and closed loop cases

% Set a PVA amplitude threshold
PVAAmpThresh = 0.1;

% Calculate the angles for the ROIs
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));

% Step through the conditions
for condID = 1:length(cond)
    
    % Make the figures
    PVADark = figure('units','normalized','outerposition',[0 0 1 1]);
    PVACL = figure('units','normalized','outerposition',[0 0 1 1]);
    
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        trialName = 'All_30C';
        
        % Step through the trials
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

            % Get the behavioral data
            tAll = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll-tAll(1);
            heading = pi/180*cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.PosRotMatch;

            % Find the dark periods
            DarkPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Direction == 0);

            % Find the closed loop period
            CLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Closed > 0);
            CLPer = setdiff(CLPer,DarkPer);
            CLPer(end) = [];            
            
            % Find periods where the stripe is visible during the closed
            % loop period
            CLVis = find(abs(heading(CLPer)) <= 2*pi/3 );
            CLVisStart = find(diff(CLPer(CLVis))>1)+1;
            CLVisStop = find(diff(CLPer(CLVis))>1);
            CLVisStart = vertcat(1,CLVisStart);
            CLVisStop = vertcat(CLVisStop,length(CLVis));
            
            % Find periods where the stripe is invisible during the closed
            % loop period
            CLInvis = find(abs(heading(CLPer)) > 2*pi/3 );
            CLInvisStart = find(diff(CLPer(CLInvis))>1)+1;
            CLInvisStop = find(diff(CLPer(CLInvis))>1);
            CLInvisStart = vertcat(1,CLInvisStart);
            CLInvisStop = vertcat(CLInvisStop,length(CLInvis));

            % Find the rotational open loop period
            OLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OLGain > 0);

            % Get the activity
            EBAct = cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;

            % Calculate the PVA
            num_ROIs = size(EBAct,1);
            clear PVA;
            clear PVAAmp;
            for ts = 1:length(EBAct)
                PVA(ts) = circ_mean(angs', squeeze(EBAct(:,ts)));
                PVAAmp(ts) = circ_r(angs', squeeze(EBAct(:,ts)));
            end

            % Keep the PVA constant if its amplitude falls below a certain
            % threshold
            PVAConst = PVA;
            for ts = 2:length(PVAConst)
                if PVAAmp(ts) < PVAAmpThresh
                    PVA(ts) = PVA(ts-1);
                end
            end
            
            % Unwrap the PVA and heading
            PVAUnwrap = UnWrap(PVA,2,0);
            headingUnwrap = UnWrap(heading,1.5,0);

            % Plot the PVA comparison in the dark
            figure(PVADark);
            subplot(cond{condID}.numFlies,3,trialID+(flyID-1)*3)
            hold on;
            plot(tAll(DarkPer),PVAUnwrap(DarkPer)-mean(PVAUnwrap(DarkPer)),'g');
            plot(tAll(DarkPer),headingUnwrap(DarkPer)-mean(headingUnwrap(DarkPer)),'k');
            xlim([tAll(DarkPer(1)) tAll(DarkPer(end))]);
            ylim([-3*pi 3*pi]);
            if flyID == 1 & trialID == 1
                text(tAll(DarkPer(1))-5,3.5*pi,cond{condID}.name,'FontSize',14);
            end
            if flyID == cond{condID}.numFlies
                xlabel('time (s)');
            end
            if trialID == 1
                ylabel('position (rad)');
            end

            % Plot the PVA comparison with a stripe
            figure(PVACL);
            subplot(cond{condID}.numFlies,3,trialID+(flyID-1)*3)
            hold on;
            plot(tAll(CLPer),PVAUnwrap(CLPer)-mean(PVAUnwrap(CLPer)),'g');
            for visSeg = 1:length(CLVisStart)
                plot(tAll(CLPer(CLVis(CLVisStart(visSeg):CLVisStop(visSeg)))),...
                    headingUnwrap(CLPer(CLVis(CLVisStart(visSeg):CLVisStop(visSeg))))-...
                    mean(headingUnwrap(CLPer)),...
                    'b');
            end
            for invisSeg = 1:length(CLInvisStart)
                plot(tAll(CLPer(CLInvis(CLInvisStart(invisSeg):CLInvisStop(invisSeg)))),...
                    headingUnwrap(CLPer(CLInvis(CLInvisStart(invisSeg):CLInvisStop(invisSeg))))-...
                    mean(headingUnwrap(CLPer)),...
                    'k');
            end
            xlim([tAll(CLPer(1)) tAll(CLPer(end))]);
            ylim([-3*pi 3*pi]);
            if flyID == 1 & trialID == 1
                text(tAll(CLPer(1))-5,3.5*pi,cond{condID}.name,'FontSize',14);
            end
            if flyID == cond{condID}.numFlies
                xlabel('time (s)');
            end
            if trialID == 1
                ylabel('position (rad)');
            end

        end
    end
end

%% Plot the PVA vs. position for the dark and closed loop cases - for dark and CL only

% Set a multiplier for the stripe
stripeMult = 360/240;

% Set a PVA amplitude threshold
PVAAmpThresh = 0;

% Calculate the angles for the ROIs
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));
    
% Step through the conditions
for condID = 1:length(cond)
    
    % Make the figures
    PVADark = figure('units','normalized','outerposition',[0 0 1 1]);
    PVACL = figure('units','normalized','outerposition',[0 0 1 1]);
    
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        trialName = 'All_30C';
        
        % Step through the trials
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

            % Get the behavioral data
            tAll = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll-tAll(1);
            heading = stripeMult*cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,2);

            % Find the dark periods
            DarkPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Direction == 0);

            % Find the closed loop period
            CLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Closed > 0);
            CLPer = setdiff(CLPer,DarkPer);
            CLPer(end) = [];            

            % Get the activity
            EBAct = cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;

            % Calculate the PVA
            num_ROIs = size(EBAct,1);
            clear PVA;
            clear PVAAmp;
            for ts = 1:length(EBAct)
                PVA(ts) = circ_mean(angs', squeeze(EBAct(:,ts)));
                PVAAmp(ts) = circ_r(angs', squeeze(EBAct(:,ts)));
            end

            % Keep the PVA constant if its amplitude falls below a certain
            % threshold
            PVAConst = PVA;
            for ts = 2:length(PVAConst)
                if PVAAmp(ts) < PVAAmpThresh
                    PVAConst(ts) = PVAConst(ts-1);
                end
            end
            
            % Clean up the PVA and heading for plotting
            PVAPlt = PVAConst;
            for tPt = 2:length(PVAPlt)
                if abs(PVAPlt(tPt)-PVAPlt(tPt-1)) > pi/2
                    PVAPlt(tPt-1) = NaN;
                end
            end
            headingPlt = heading;
            for tPt = 2:length(headingPlt)
                if abs(headingPlt(tPt)-headingPlt(tPt-1)) > pi
                    headingPlt(tPt-1) = NaN;
                end
            end
                    
            
            % Plot the PVA comparison in the dark
            figure(PVADark);
            subplot(6,1,trialID+(flyID-1)*3)
            hold on;
            imagesc(tAll(DarkPer),angs,EBAct(:,DarkPer));
            colormap(brewermap(64, 'Blues'));
            xlim([tAll(DarkPer(1)) tAll(DarkPer(end))]);
            ylim([-pi pi]);
%             plot(tAll(DarkPer),PVAPlt(DarkPer),'Color',[0 0 0.5]);
            plot(tAll(DarkPer),headingPlt(DarkPer),'k');
            if flyID == 1 & trialID == 1
                text(tAll(DarkPer(1))-5,1.75*pi,strcat(cond{condID}.name,'-dark'),'FontSize',14);
            end
            if flyID == cond{condID}.numFlies & trialID == 3
                xlabel('time (s)');
                ylabel('position (rad)');
            end
            caxis([0 1]);
            colorbar;
            text(tAll(DarkPer(1))-2,1.25*pi,strcat('fly #',num2str(flyID),'-trial #',num2str(trialID)));

            % Plot the PVA comparison with a stripe
            figure(PVACL);
            subplot(6,1,trialID+(flyID-1)*3)
            hold on;
            imagesc(tAll(CLPer),angs,EBAct(:,CLPer));
            colormap(brewermap(64, 'Blues'));
%             plot(tAll(CLPer),PVAPlt(CLPer),'Color',[0 0 0.5]);
            plot(tAll(CLPer),...
                headingPlt(CLPer),...
                'b');
            xlim([tAll(CLPer(1)) tAll(CLPer(end))]);
            ylim([-pi pi]);
            if flyID == 1 & trialID == 1
                text(tAll(CLPer(1))-5,1.75*pi,strcat(cond{condID}.name,'-CL'),'FontSize',14);
            end
            if flyID == cond{condID}.numFlies & trialID == 3
                xlabel('time (s)');
                ylabel('position (rad)');
            end
            caxis([0 1]);
            colorbar;
            text(tAll(CLPer(1))-4,1.25*pi,strcat('fly #',num2str(flyID),'-trial #',num2str(trialID)));

        end
    end
    
%     set(PVADark,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%     print(PVADark,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\Shi\Figures\ShiExample_',cond{condID}.name,'_Dark'),'-dpdf');
%     set(PVACL,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%     print(PVACL,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\Shi\Figures\ShiExample_',cond{condID}.name,'_CL'),'-dpdf');
end

%% Filter out the high frequency components and compare the heading and PVA

% Set a PVA amplitude threshold
PVAAmpThresh = 0.2;

% Calculate the angles for the ROIs
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));

% Build a Butterworth filter for smoothing the data
order = 6; % Order of the filter
lco = 0.0001; % Low pass frequency
hco = 0.5; % High pass frequency
Fs = 8.5; % Sampling frequency\
n=round(order/2);
wn=[lco hco]/(Fs/2);
[b,a]=butter(n, wn);


% Step through the conditions
for condID = 1:length(cond)
    
    % Make the figures
    PVADark = figure('units','normalized','outerposition',[0 0 1 1]);
    PVACL = figure('units','normalized','outerposition',[0 0 1 1]);
    
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        trialName = 'All_30C';
        
        % Step through the trials
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

            % Get the behavioral data
            tAll = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll-tAll(1);
            heading = pi/180*cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.PosRotMatch;

            % Find the dark periods
            DarkPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Direction == 0);

            % Find the closed loop period
            CLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Closed > 0);
            CLPer = setdiff(CLPer,DarkPer);
            CLPer(end) = [];            
            
            % Find the rotational open loop period
            OLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OLGain > 0);

            % Get the activity
            EBAct = cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;

            % Calculate the PVA
            num_ROIs = size(EBAct,1);
            clear PVA;
            clear PVAAmp;
            for ts = 1:length(EBAct)
                PVA(ts) = circ_mean(angs', squeeze(EBAct(:,ts)));
                PVAAmp(ts) = circ_r(angs', squeeze(EBAct(:,ts)));
            end

            % Keep the PVA constant if its amplitude falls below a certain
            % threshold
            PVAConst = PVA;
            for ts = 2:length(PVAConst)
                if PVAAmp(ts) < PVAAmpThresh
                    PVA(ts) = PVA(ts-1);
                end
            end
            
            % Unwrap the PVA and heading
            PVAUnwrap = UnWrap(PVA,1.5,0);
            headingUnwrap = UnWrap(heading,1.5,0);

            % Plot the PVA comparison in the dark
            figure(PVADark);
            subplot(cond{condID}.numFlies,3,trialID+(flyID-1)*3)
            hold on;
%             plot(tAll(DarkPer),PVAUnwrap(DarkPer)-mean(PVAUnwrap(DarkPer)),'g');
            plot(tAll(DarkPer),filtfilt(b,a,PVAUnwrap(DarkPer))-mean(filtfilt(b,a,PVAUnwrap(DarkPer))),'g');
%             plot(tAll(DarkPer),headingUnwrap(DarkPer)-mean(headingUnwrap(DarkPer)),'k');
            plot(tAll(DarkPer),filtfilt(b,a,headingUnwrap(DarkPer))-mean(filtfilt(b,a,headingUnwrap(DarkPer))),'k');
            xlim([tAll(DarkPer(1)) tAll(DarkPer(end))]);
            ylim([-3*pi 3*pi]);
            if flyID == 1 & trialID == 1
                text(tAll(DarkPer(1))-5,3.5*pi,cond{condID}.name,'FontSize',14);
            end
            if flyID == cond{condID}.numFlies
                xlabel('time (s)');
            end
            if trialID == 1
                ylabel('position (rad)');
            end

            % Plot the PVA comparison with a stripe
            figure(PVACL);
            subplot(cond{condID}.numFlies,3,trialID+(flyID-1)*3)
            hold on;
%             plot(tAll(CLPer),PVAUnwrap(CLPer)-mean(PVAUnwrap(CLPer)),'g');
            plot(tAll(CLPer),filtfilt(b,a,PVAUnwrap(CLPer))-mean(filtfilt(b,a,PVAUnwrap(CLPer))),'g');
%             plot(tAll(CLPer),headingUnwrap(CLPer)-mean(headingUnwrap(CLPer)),'b');
            plot(tAll(CLPer),filtfilt(b,a,headingUnwrap(CLPer))-mean(filtfilt(b,a,headingUnwrap(CLPer))),'b');
            xlim([tAll(CLPer(1)) tAll(CLPer(end))]);
            ylim([-3*pi 3*pi]);
            if flyID == 1 & trialID == 1
                text(tAll(CLPer(1))-5,3.5*pi,cond{condID}.name,'FontSize',14);
            end
            if flyID == cond{condID}.numFlies
                xlabel('time (s)');
            end
            if trialID == 1
                ylabel('position (rad)');
            end

        end
    end
end

%% Quantify the above - look at as a function of SD
% Specify correlation and S.D. bins
corrEdges = linspace(-1,1,21);
corrCents = corrEdges(1:end-1)+0.5*mean(diff(corrEdges));
errEdges = linspace(0,100,16);
errCents = errEdges(1:end-1)+0.5*mean(diff(errEdges));
SDEdges = linspace(0,2.5,18);
SDCents = SDEdges(1:end-1)+0.5*mean(diff(SDEdges));

% Specify the sliding window parameters
frameNum = 32;
frameStep = 8;

% Set a PVA amplitude threshold
PVAAmpThresh = 0.2;

% Calculate the angles for the ROIs
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));

% Chunk the data into windows
% Step through the conditions
for condID = 1:length(cond)
    
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        % Step through the periods
        for periodID = 1:4
        
            % Step through the trials
            for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

                % Get the behavioral data
                heading = pi/180*allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.heading;
                
                % Get the activity
                EBAct = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.act-1;
                
                % Calculate the PVA
                num_ROIs = size(EBAct,1);
                clear PVA;
                clear PVAAmp;
                for ts = 1:length(EBAct)
                    PVA(ts) = circ_mean(angs', squeeze(EBAct(:,ts)));
                    PVAAmp(ts) = circ_r(angs', squeeze(EBAct(:,ts)));
                end

                % Keep the PVA constant if its amplitude falls below a certain
                % threshold
                PVAConst = PVA;
                for ts = 2:length(PVAConst)
                    if PVAAmp(ts) < PVAAmpThresh
                        PVA(ts) = PVA(ts-1);
                    end
                end

                % Unwrap the PVA and heading
                PVAUnwrap = UnWrap(PVA,1.5,0);
                headingUnwrap = UnWrap(heading,1.5,0);
                
                % Chunk the data into windows
                numWins = floor((length(PVAUnwrap)-frameNum)/frameStep);
                winData = zeros(2,numWins,frameNum);
                for winNow = 1:numWins
                    winData(1,winNow,:) = PVAUnwrap(1+frameStep*(winNow-1):frameNum+frameStep*(winNow-1))-PVAUnwrap(1+frameStep*(winNow-1));
                    winData(2,winNow,:) = headingUnwrap(1+frameStep*(winNow-1):frameNum+frameStep*(winNow-1))-headingUnwrap(1+frameStep*(winNow-1));
                end
                allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.winData = winData;
            end
        end
    end
end

% Find the heading distribution
headingAll = [];
for condID = 1:length(cond)
    
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        % Step through the periods
        for periodID = 1:4
        
            % Step through the trials
            for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

                heading = squeeze(allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.winData(2,:,:));
                for binNow = 1:size(heading,1)
                    headingAll = horzcat(headingAll,heading(binNow,:));
                end
            end
        end
    end
end
headingSD = std(headingAll(find(headingAll~=0)));

% Find the correlations and bin the data accordingly

matCorr = zeros(length(cond),4,length(corrEdges)-1,length(SDEdges)-1);
matErr = zeros(length(cond),4,length(corrEdges)-1,length(SDEdges)-1);
for condID = 1:length(cond)
    
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        % Step through the periods
        for periodID = 1:4
            winData = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.winData;
            for binNow = 1:size(winData,2)
                PVANow = winData(1,binNow,:);
                headingNow = winData(2,binNow,:);
                errorNow = sum(abs(PVANow-headingNow));
                CCs = corrcoef(PVANow,headingNow);
                SDBin = discretize(std(headingNow)/headingSD,SDEdges);
                errBin = discretize(errorNow,errEdges);
                if ~isnan(errBin)
                    matErr(condID,periodID,errBin,SDBin) = matErr(condID,periodID,errBin,SDBin) + 1;
                end
                if ~isnan(CCs(2,1))
                    CCBin = discretize(CCs(2,1),corrEdges);
                    matCorr(condID,periodID,CCBin,SDBin) = matCorr(condID,periodID,CCBin,SDBin) + 1;
                end
            end
        end
    end
end

% Plot the matrices
figure;
for condID = 1:length(cond)
    for periodID = 1:4
        subplot(length(cond),4,periodID+4*(condID-1));
        imagesc(corrCents,SDCents,squeeze(matCorr(condID,periodID,:,:))');
    end
end

% Plot the matrices
figure;
for condID = 1:length(cond)
    for periodID = 1:4
        subplot(length(cond),4,periodID+4*(condID-1));
        imagesc(errCents,SDCents,squeeze(matErr(condID,periodID,:,:))');
    end
end

%% Look at the bump amplitude vs. forward or rotational velocity as a function of condition

% Specify bins for velocities
vRBinEdges = linspace(0,4*pi,21);
vFBinEdges = linspace(-0.1,2,11);

% Make the figures
ampComp = figure('units','normalized','outerposition',[0 0 1 1]);

colorCond = zeros(3,3);
colorCond(1,:) = [0 0 0];
colorCond(2,:) = [0.25 0.25 1];
colorCond(3,:) = [0.5 0.25 0];
    
% Step through the conditions
for condID = [1 3];%1:length(cond)
    
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        trialName = 'All_30C';
        
        maxAllDark = [];
        maxAllCL = [];
        maxAllOL = [];
        
        % Step through the trials
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

            % Get the behavioral data
            tAll = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll-tAll(1);
            vR = abs(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vRot);
            vF = abs(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vF);

            % Find the dark periods
            DarkPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Direction == 0);

            % Find the closed loop period
            CLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Closed > 0);
            CLPer = setdiff(CLPer,DarkPer);
            CLPer(end) = [];            
            
            % Find the rotational open loop period
            OLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OLGain > 0);
            OLPer(end) = [];

            % Get the activity max
            EBActMax = max(cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1);

            maxAllDark = horzcat(maxAllDark,...
                vertcat(vertcat(vR(DarkPer)',vF(DarkPer)'),EBActMax(DarkPer)));
            maxAllCL = horzcat(maxAllCL,...
                vertcat(vertcat(vR(CLPer)',vF(CLPer)'),EBActMax(CLPer)));
            maxAllOL = horzcat(maxAllOL,...
                vertcat(vertcat(vR(OLPer)',vF(OLPer)'),EBActMax(OLPer)));

        end
        
        vRIDsDark = discretize(maxAllDark(1,:),vRBinEdges);
        vFIDsDark = discretize(maxAllDark(2,:),vFBinEdges);
        vRIDsCL = discretize(maxAllCL(1,:),vRBinEdges);
        vFIDsCL = discretize(maxAllCL(2,:),vFBinEdges);
        vRIDsOL = discretize(maxAllOL(1,:),vRBinEdges);
        vFIDsOL = discretize(maxAllOL(2,:),vFBinEdges);
        
        for stim = 1:3
            if stim == 1
                vRNow = vRIDsDark;
                vFNow = vFIDsDark;
                maxNow = maxAllDark(3,:);
                titleNow = 'dark';
            elseif stim == 2
                vRNow = vRIDsCL;
                vFNow = vFIDsCL;
                maxNow = maxAllCL(3,:);
                titleNow = 'closed loop';
            else
                vRNow = vRIDsOL;
                vFNow = vFIDsOL;
                maxNow = maxAllOL(3,:);
                titleNow = 'open loop';
            end
            subplot(2,3,stim);
            hold on;
            for binID = 1:length(vRBinEdges)-1
                ptsNow = find(vRNow == binID);
                if ~isempty(ptsNow)
                    patch('XData',[vRBinEdges(binID) vRBinEdges(binID+1) vRBinEdges(binID+1) vRBinEdges(binID)],...
                        'YData',[0 0 mean(maxNow(ptsNow)) mean(maxNow(ptsNow))],...
                        'EdgeColor',colorCond(condID,:),'FaceColor',colorCond(condID,:),'FaceAlpha',0.2);
    %                 rectangle('Position',[vRBinEdges(binID) 0 vRBinEdges(binID+1)-vRBinEdges(binID) ...
    %                     mean(maxAllDark(3,ptsNow))],...
    %                     'EdgeColor',colorCond(condID,:));
                end
            end
            xlim([0 2*pi]);
            xlabel('vR (rad/s)');
            ylim([0 2]);
            if stim == 1
                ylabel('max DF/F');
                if flyID == 1
                    patch('XData',[1.75*pi-0.1 1.75*pi+0.1 1.75*pi+0.1 1.75*pi-0.1],...
                        'YData',[2-0.1*condID 2-0.1*condID 2-0.1*(condID-1) 2-0.1*(condID-1)],...
                        'EdgeColor',colorCond(condID,:),'FaceColor',colorCond(condID,:),'FaceAlpha',0.2);
                    text(1.75*pi+0.15,2-0.1*(condID-0.5),strcat(cond{condID}.name,' (n = ',num2str(cond{condID}.numFlies),')'));
                end
            end
            title(titleNow);
            
            
            subplot(2,3,stim+3);
            hold on;
            for binID = 1:length(vFBinEdges)-1
                ptsNow = find(vFNow == binID);
                if ~isempty(ptsNow)
                    patch('XData',[vFBinEdges(binID) vFBinEdges(binID+1) vFBinEdges(binID+1) vFBinEdges(binID)],...
                        'YData',[0 0 mean(maxNow(ptsNow)) mean(maxNow(ptsNow))],...
                        'EdgeColor',colorCond(condID,:),'FaceColor',colorCond(condID,:),'FaceAlpha',0.2);
    %                 rectangle('Position',[vRBinEdges(binID) 0 vRBinEdges(binID+1)-vRBinEdges(binID) ...
    %                     mean(maxAllDark(3,ptsNow))],...
    %                     'EdgeColor',colorCond(condID,:));
                end
            end
            xlim([0 2]);
            xlabel('vF (cm/s)');
            ylim([0 2]);
            if stim == 1
                ylabel('max DF/F');
            end
            
        end
        
        
    end
end

% set(ampComp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(ampComp,'C:\Users\turnerevansd\Documents\RawAnalysis\Shi\Figures\BumpAmplitudes_GEVsCtrl','-dpdf');

%% Look at the bump amplitude vs. forward or rotational velocity as a function of condition - V2

% Savitsky-Golay filter parameters
sgoOrder = 3;
sgoFrames = 11;

% Specify bins for velocities
vRBinEdges = linspace(0,4*pi,41);
vFBinEdges = linspace(-0.2,2,23);

% Make the figures
ampComp = figure('units','normalized','outerposition',[0 0 1 1]);

colorCond = zeros(3,3);
colorCond(1,:) = [0 0 0];
colorCond(2,:) = [0.25 0.25 1];
colorCond(3,:) = [0.5 0.25 0];
    
% Step through the conditions
for condID = 1:length(cond)
    
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        trialName = 'All_30C';
        
        maxAllDark = [];
        maxAllCL = [];
        maxAllOL = [];
        
        % Step through the trials
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))

            % Get the behavioral data
            tAll = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll-tAll(1);
            vR = abs(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vRot);
            vR = sgolayfilt(vR,sgoOrder,sgoFrames);
            vF = abs(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vF);
            vF = sgolayfilt(vF,sgoOrder,sgoFrames);
            
            % Find the dark periods
            DarkPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Direction == 0);

            % Find the closed loop period
            CLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Closed > 0);
            CLPer = setdiff(CLPer,DarkPer);
            CLPer(end) = [];            
            
            % Find the rotational open loop period
            OLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OLGain > 0);
            OLPer(end) = [];

            % Get the activity max
            EBActMax = max(cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1);

            maxAllDark = horzcat(maxAllDark,...
                vertcat(vertcat(vR(DarkPer)',vF(DarkPer)'),EBActMax(DarkPer)));
            maxAllCL = horzcat(maxAllCL,...
                vertcat(vertcat(vR(CLPer)',vF(CLPer)'),EBActMax(CLPer)));
            maxAllOL = horzcat(maxAllOL,...
                vertcat(vertcat(vR(OLPer)',vF(OLPer)'),EBActMax(OLPer)));

        end
        
        vRIDsDark = discretize(maxAllDark(1,:),vRBinEdges);
        vFIDsDark = discretize(maxAllDark(2,:),vFBinEdges);
        vRIDsCL = discretize(maxAllCL(1,:),vRBinEdges);
        vFIDsCL = discretize(maxAllCL(2,:),vFBinEdges);
        vRIDsOL = discretize(maxAllOL(1,:),vRBinEdges);
        vFIDsOL = discretize(maxAllOL(2,:),vFBinEdges);
        
        for stim = 1:3
            if stim == 1
                vRNow = maxAllDark(1,:);
                vRDisc = vRIDsDark;
                vFNow = maxAllDark(2,:);
                vFDisc = vFIDsDark;
                maxNow = maxAllDark(3,:);
                titleNow = 'dark';
            elseif stim == 2
                vRNow = maxAllCL(1,:);
                vRDisc = vRIDsCL;
                vFNow = maxAllCL(2,:);
                vFDisc = vFIDsCL;
                maxNow = maxAllCL(3,:);
                titleNow = 'closed loop';
            else
                vRNow = maxAllOL(1,:);
                vRDisc = vRIDsOL;
                vFNow = maxAllOL(2,:);
                vFDisc = vFIDsOL;
                maxNow = maxAllOL(3,:);
                titleNow = 'open loop';
            end
            subplot(2,3,stim);
            hold on;
%             scatter(vRNow,maxNow,20,colorCond(condID,:),'filled');
            for binID = 1:length(vRBinEdges)-1
                ptsNow = find(vRDisc == binID);
                if ~isempty(ptsNow)
                    binCent = 0.5*(vRBinEdges(binID+1)+vRBinEdges(binID))+0.02*flyID+0.1*(condID-1);
%                     scatter(binCent,mean(maxNow(ptsNow)),40,colorCond(condID,:),'filled');
%                     line([binCent binCent],...
%                         [mean(maxNow(ptsNow))-std(maxNow(ptsNow)) mean(maxNow(ptsNow))+std(maxNow(ptsNow))],...
%                         'Color',colorCond(condID,:));
                    patch('XData',[binCent-0.02 binCent+0.02 binCent+0.02 binCent-0.02],...
                        'YData',[mean(maxNow(ptsNow))-std(maxNow(ptsNow)) mean(maxNow(ptsNow))-std(maxNow(ptsNow)) ...
                        mean(maxNow(ptsNow))+std(maxNow(ptsNow)) mean(maxNow(ptsNow))+std(maxNow(ptsNow))],...
                        'EdgeColor',colorCond(condID,:),'FaceColor',colorCond(condID,:),'FaceAlpha',0.2);
                end
            end
%             alpha(0.2);
            xlim([0 pi/2+0.15]);
            xlabel('vR (rad/s)');
            ylim([0 2]);
            if stim == 1
                ylabel('max DF/F');
                if flyID == 1
                    patch('XData',[0.5*pi-0.05 0.5*pi+0.05 0.5*pi+0.05 0.5*pi-0.05],...
                        'YData',[2-0.15*condID 2-0.15*condID 1.9-0.15*(condID-1) 1.9-0.15*(condID-1)],...
                        'EdgeColor',colorCond(condID,:),'FaceColor',colorCond(condID,:),'FaceAlpha',0.2);
                    text(0.5*pi+0.15,1.95-0.15*(condID-0.5),strcat(cond{condID}.name,' (n = ',num2str(cond{condID}.numFlies),')'));
                end
            end
            title(titleNow);
            
            
            subplot(2,3,stim+3);
            hold on;
%             scatter(vFNow,maxNow,20,colorCond(condID,:),'filled');
            for binID = 1:length(vFBinEdges)-1
                ptsNow = find(vFDisc == binID);
                if ~isempty(ptsNow)
                    binCent = 0.5*(vFBinEdges(binID+1)+vFBinEdges(binID))+0.005*flyID+0.025*(condID-2);
%                     scatter(binCent,mean(maxNow(ptsNow)),40,colorCond(condID,:),'filled');
%                     line([binCent binCent],...
%                         [mean(maxNow(ptsNow))-std(maxNow(ptsNow)) mean(maxNow(ptsNow))+std(maxNow(ptsNow))],...
%                         'Color',colorCond(condID,:));
                    patch('XData',[binCent-0.005 binCent+0.005 binCent+0.005 binCent-0.005],...
                        'YData',[mean(maxNow(ptsNow))-std(maxNow(ptsNow)) mean(maxNow(ptsNow))-std(maxNow(ptsNow)) ...
                        mean(maxNow(ptsNow))+std(maxNow(ptsNow)) mean(maxNow(ptsNow))+std(maxNow(ptsNow))],...
                        'EdgeColor',colorCond(condID,:),'FaceColor',colorCond(condID,:),'FaceAlpha',0.2);
                end
            end
%             alpha(0.2);
            xlim([0 0.6]);
            xlabel('vF (cm/s)');
            ylim([0 2]);
            if stim == 1
                ylabel('max DF/F');
            end
            
        end
        
        
    end
end

set(ampComp,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(ampComp,'C:\Users\turnerevansd\Documents\RawAnalysis\Shi\Figures\BumpAmplitudes_All_RT','-dpdf');

%% Look at the bump persistance
condID = 1;

% Set the Savitsky Golay filter parameters
sgolayOrder = 3;
sgolayFrames = 11;

% Set the velocity thresholds
vRThresh = 0.01*pi;
vFThresh = 0.01;

% Set a threshold for the length of a standing bout
minBoutL = 5;

% Step through the flies
for flyID = 1:cond{condID}.numFlies

    % Sort the data by period
    for periodID = 1:4
        
        numPeriods = length(allAct{condID}.fly{flyID}.period);
        maxBoutL = 0; % Max bout length
        numBouts = 0; % Number of bouts
            
        % Sort the data across trials
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.All_30C)

            % Find the time step
            tPts = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.tPts;
            tStep = mean(diff(tPts));

            % Find when the animal is standing
            vR = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vR;
            vR = sgolayfilt(vR,sgolayOrder,sgolayFrames);
            vF = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vF;
            vF = sgolayfilt(vF,sgolayOrder,sgolayFrames);
            standingBouts = intersect(find(vR<=vRThresh),find(vF<=vFThresh));

            if ~isempty(standingBouts)
                % Group the standing bouts
                boutStartIndex = vertcat(1,find(diff(standingBouts)>1)+1);
                boutStopIndex = vertcat(find(diff(standingBouts)>1),length(standingBouts));

                % Remove shorter bouts
                shortBouts = [];
                for bout = 1:length(boutStartIndex)
                 if (standingBouts(boutStopIndex(bout))-standingBouts(boutStartIndex(bout))) <= minBoutL 
                     shortBouts = vertcat(shortBouts,bout);
                 end
                end
                boutStartIndex(shortBouts) = [];
                boutStopIndex(shortBouts) = [];

                maxBoutL = max(maxBoutL,max(boutStopIndex-boutStartIndex));
                numBouts = numBouts+length(boutStartIndex);

                % Add the standing bout info to the object
                allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.boutStart = standingBouts(boutStartIndex);
                allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.boutStop = standingBouts(boutStopIndex);
            end

        end
        
        allAct{condID}.fly{flyID}.period{periodID}.maxBoutL = maxBoutL;
        allAct{condID}.fly{flyID}.period{periodID}.numBouts = numBouts;
        
    end
end

% Group the standing bouts across trials

% Make the figure
bumpPer = figure('units','normalized','outerposition',[0 0 1 1]);

% Step through the flies
for flyID = 1:cond{condID}.numFlies
    
    % Sort the data by period
    for periodID = 1:4
    
        numPeriods = length(allAct{condID}.fly{flyID}.period);
        maxBoutL = allAct{condID}.fly{flyID}.period{periodID}.maxBoutL;
        numBouts = allAct{condID}.fly{flyID}.period{periodID}.numBouts;
        
        standAct = zeros(numBouts,16,maxBoutL);
        standActMean = zeros(16,maxBoutL);
        standBout = 1;

        % Sort the data across trials
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.All_30C)

            % Place all of the standing bout activity into the array
            actNow = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.act-1;
            if ~isfield(allAct{condID}.fly{flyID}.period{periodID}.trial{trialID},'boutStart');
                continue;
            end
                
            boutStart = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.boutStart;
            boutStop = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.boutStop;
            
            if ~isempty(boutStart)
                for bout = 1:length(boutStart)
                    actBout = actNow(:,boutStart(bout):boutStop(bout));
                    maxROI = find(actBout(:,1) == max(actBout(:,1)));
                    actBoutShift = circshift(actBout,8-maxROI,1);
                    standAct(standBout,:,1:(boutStop(bout)-boutStart(bout)+1)) = actBoutShift;
                    standBout = standBout + 1;
                end
            end
        end
        
        for l = 1:maxBoutL
            allActNow = squeeze(standAct(:,:,l));
            zeroPers = [];
            for bout = 1:size(allActNow,1)
                if max(allActNow(bout,:)) == 0
                    zeroPers = vertcat(zeroPers,bout);
                end
            end
            allActNow(zeroPers,:) = [];
            standActMean(:,l) = mean(allActNow,1);
        end
        
        subplot(numPeriods,cond{condID}.numFlies,flyID+cond{condID}.numFlies*(periodID-1))
        tVals = tStep*[0:size(standActMean,2)-1];
        angVals = linspace(-pi,pi,17);
        angVals(end) = [];
        angVals = angVals + 0.5*mean(diff(angVals));
        imagesc(tVals,angVals,standActMean);
        xlim([0 6]);
        caxis([0 0.75]);
        colormap(brewermap(64, 'Blues'));
        if periodID == 4
            xlabel('time (s)');
        end
        if flyID == 1
            ylabel('EB position (rad');
        end
        if periodID == 1
           text(0,-1.25*pi,strcat('fly #',num2str(flyID))); 
        end
        title(allAct{condID}.fly{flyID}.period{periodID}.type);
    end
end

% set(bumpPer,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(bumpPer,'C:\Users\turnerevansd\Documents\RawAnalysis\G-E\BumpPersistance','-dpdf');

%% Look at bump shape during turns

% Specify the bins
vRBins = linspace(-1.5*pi,1.5*pi,7);
vRCents = vRBins(1:end-1)+0.5*mean(diff(vRBins));

% Savitsky-Golay filter the rotational velocities
sgoOrder = 3;
sgoFrames = 11;

% Group the cross-sections by velocity
% Step through the conditions
for condID = 1:length(cond)
    
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        % Step through the periods
        for periodID = 1:4
        
            % Step through the trials
            for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))
                
                % Get the activity and rotational velocities
                act = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.act;
                vR = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vR;
                vR = sgolayfilt(vR,sgoOrder,sgoFrames);
                
                % Circularly shift the activity
                for tPt = 1:length(act)
                    actNow = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.act(:,tPt);
                    maxAct = max(actNow);
                    pkPos = find(actNow == maxAct);
                    act(:,tPt) = circshift(actNow,8-pkPos);
                end
                
                % Discretize the rotational velocity
                vRDisc = discretize(vR,vRBins);
                
                % Group the activity profiles by velocity
                if trialID == 1
                    allAct{condID}.fly{flyID}.period{periodID}.actSort = cell(length(vRCents),1);
                end
                for vRStep = 1:length(vRCents)
                    actSort = act(:,find(vRDisc == vRStep));
                    allAct{condID}.fly{flyID}.period{periodID}.actSort{vRStep} = ...
                        horzcat(allAct{condID}.fly{flyID}.period{periodID}.actSort{vRStep},actSort);
                end
            end
        end
    end
end

% Plot the mean and SD profiles for each condition, each fly, and each
% stimulus
turnColors = jet(length(vRCents));
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));
% Step through the conditions
for condID = 1:length(cond)
    
    profPlot = figure('units','normalized','outerposition',[0 0 1 1]);
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        % Step through the periods
        for periodID = 1:4
            subplot(5,4,periodID+4*(flyID-1));
            hold on;
            for vRStep = [1 2 5 6]
                if ~isempty(allAct{condID}.fly{flyID}.period{periodID}.actSort{vRStep})
                    plot(angs,...
                        mean(allAct{condID}.fly{flyID}.period{periodID}.actSort{vRStep},2),...
                        'Color',turnColors(vRStep,:));
                    ylim([1 2]);
                end
            end
        end
    end
end



%% Plot all dark and CL periods

condID = 1;

PVADark = figure('units','normalized','outerposition',[0 0 1 1]);
PVACL = figure('units','normalized','outerposition',[0 0 1 1]);

angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));

for flyID = 1:cond{condID}.numFlies
    for hc = 1:2
        if hc == 1
            trialName = 'All';
        else
            trialName = 'All_30C';
        end
        meanIntAllDark = [];
        meanIntAllCL = [];
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName));
            
            % Get the behavioral data
            tAll = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll-tAll(1);
            heading = pi/180*cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.PosRotMatch;
            
            % Find the dark periods
            DarkPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Direction == 0);

            % Find the closed loop period
            CLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Closed > 0);
            CLPer = setdiff(CLPer,DarkPer);
            CLPer(end) = [];

            % Find the rotational open loop period
            OLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OLGain > 0);
            
            % Get the activity
            EBAct = cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;

            % Plot the PVA comparison in the dark
            figure(PVADark);
            subplot(2*cond{condID}.numFlies,3,trialID+(hc-1)*3+(flyID-1)*6)
            hold on;
            imagesc(tAll(DarkPer),angs,EBAct(:,DarkPer));
            plot(tAll(DarkPer),heading(DarkPer),'k');
            xlim([tAll(DarkPer(1)) tAll(DarkPer(end))]);
            ylim([-pi pi]);
            colormap(brewermap(64, 'Blues'));
            caxis([0 2]);
            
            % Plot the PVA comparison with a stripe
            figure(PVACL);
            subplot(2*cond{condID}.numFlies,3,trialID+(hc-1)*3+(flyID-1)*6)
            hold on;
            imagesc(tAll(CLPer),angs,EBAct(:,CLPer));
            plot(tAll(CLPer),heading(CLPer),'b');
            xlim([tAll(CLPer(1)) tAll(CLPer(end))]);
            ylim([-pi pi]); 
            colormap(brewermap(64, 'Blues'));
            caxis([0 2]);
        end
    end
end

%% Plot the PVA vs. position for the dark and closed loop cases

condID =3;

PVADark = figure('units','normalized','outerposition',[0 0 1 1]);
PVACL = figure('units','normalized','outerposition',[0 0 1 1]);

angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));

meanIntBins = linspace(0,1,21);
meanIntCents = meanIntBins+0.5*mean(diff(meanIntBins));
meanIntCents(end) = [];

for flyID = 1:cond{condID}.numFlies
    for hc = 1:2
        if hc == 1
            trialName = 'All';
        else
            trialName = 'All_30C';
        end
        meanIntAllDark = [];
        meanIntAllCL = [];
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))
            
            % Get the behavioral data
            tAll = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll-tAll(1);
            heading = pi/180*cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.PosRotMatch;
            
            % Find the dark periods
            DarkPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Direction == 0);

            % Find the closed loop period
            CLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Closed > 0);
            CLPer = setdiff(CLPer,DarkPer);
            CLPer(end) = [];

            % Find the rotational open loop period
            OLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OLGain > 0);
            
            % Get the activity
            EBAct = cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;
            
            % Calculate the PVA
            num_ROIs = size(EBAct,1);
            clear meanAngRaw;
            clear meanIntRaw;
            for ts = 1:length(EBAct)
                meanAngRaw(ts) = circ_mean(angs', squeeze(EBAct(:,ts)));
                meanIntRaw(ts) = circ_r(angs', squeeze(EBAct(:,ts)));
            end
            meanIntAllDark = horzcat(meanIntAllDark,meanIntRaw(DarkPer));
            meanIntAllCL = horzcat(meanIntAllCL,meanIntRaw(CLPer));
            
            PVAPlt = meanAngRaw;
            for ts = 1:length(PVAPlt)-1
                if abs(PVAPlt(ts+1) - PVAPlt(ts)) > pi
                    PVAPlt(ts) = NaN; 
                end
            end

            % Plot the PVA comparison in the dark
            figure(PVADark);
            subplot(2*cond{condID}.numFlies,4,trialID+(hc-1)*4+(flyID-1)*8)
            hold on;
            plot(tAll(DarkPer),PVAPlt(DarkPer),'g');
            plot(tAll(DarkPer),heading(DarkPer),'k');
            xlim([tAll(DarkPer(1)) tAll(DarkPer(end))]);
            ylim([-pi pi]);
            
            % Plot the PVA comparison with a stripe
            figure(PVACL);
            subplot(2*cond{condID}.numFlies,4,trialID+(hc-1)*4+(flyID-1)*8)
            hold on;
            plot(tAll(CLPer),PVAPlt(CLPer),'g');
            plot(tAll(CLPer),heading(CLPer),'b');
            xlim([tAll(CLPer(1)) tAll(CLPer(end))]);
            ylim([-pi pi]);
            
        end
        figure(PVADark);
        subplot(2*cond{condID}.numFlies,4,4+(hc-1)*4+(flyID-1)*8)
        histogram(meanIntAllDark,meanIntBins);
        xlim([0 1]);
        ylim([0 300]);
        
        figure(PVACL);
        subplot(2*cond{condID}.numFlies,4,4+(hc-1)*4+(flyID-1)*8)
        histogram(meanIntAllCL,meanIntBins);
        xlim([0 1]);
        ylim([0 300]);
    end
end

%% Plot the bump max vs. vR and vF

condID = 2;

maxDark = figure('units','normalized','outerposition',[0 0 1 1]);
maxCL = figure('units','normalized','outerposition',[0 0 1 1]);

angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));

for flyID = 1:cond{condID}.numFlies
    for hc = 1:2
        if hc == 1
            trialName = 'All';
        else
            trialName = 'All_30C';
        end
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName));
            
            % Get the behavioral data
            tAll = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll-tAll(1);
            vR = abs(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vRot);
            vF = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.vF;
            
            % Find the dark periods
            DarkPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Direction == 0);

            % Find the closed loop period
            CLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Closed > 0);
            CLPer = setdiff(CLPer,DarkPer);
            CLPer(end) = [];

            % Find the rotational open loop period
            OLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OLGain > 0);
            
            % Get the activity
            EBActMax = max(cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1);
            EBActMax(end) = [];

            % Plot the PVA comparison in the dark
            figure(maxDark);
            subplot(2*cond{condID}.numFlies,2,2*flyID-1)
            hold on;
            if hc == 1
                scatter(vR(DarkPer),EBActMax(DarkPer),50,'b','filled');
            else
                scatter(vR(DarkPer),EBActMax(DarkPer),50,'r','filled');
            end
            alpha(0.2);
            xlabel('vR (rad/s)');
            xlim([0 2*pi]);
            ylim([0 2.5]);
            
            subplot(2*cond{condID}.numFlies,2,2*flyID)
            hold on;
            if hc == 1
                scatter(vF(DarkPer),EBActMax(DarkPer),50,'b','filled');
            else
                scatter(vF(DarkPer),EBActMax(DarkPer),50,'r','filled');
            end
            alpha(0.2);
            xlabel('vF (cm/s)');
            xlim([0 2]);
            ylim([0 2.5]);
            
            % Plot the PVA comparison with a stripe
            figure(maxCL);
            subplot(2*cond{condID}.numFlies,2,2*flyID-1)
            hold on;
            if hc == 1
                scatter(vR(CLPer),EBActMax(CLPer),50,'b','filled');
            else
                scatter(vR(CLPer),EBActMax(CLPer),50,'r','filled');
            end
            alpha(0.2);
            xlabel('vR (rad/s)');
            xlim([0 2*pi]);
            ylim([0 2.5]);
            
            subplot(2*cond{condID}.numFlies,2,2*flyID)
            hold on;
            if hc == 1
                scatter(vF(CLPer),EBActMax(CLPer),50,'b','filled');
            else
                scatter(vF(CLPer),EBActMax(CLPer),50,'r','filled');
            end
            alpha(0.2);
            xlabel('vF (cm/s)');
            xlim([0 2]);
            ylim([0 2.5]);
            
        end
    end
end

%% Plot the PVA vs. position for the closed loop, high temperature case

PVADark = figure('units','normalized','outerposition',[0 0 1 1]);
PVACL = figure('units','normalized','outerposition',[0 0 1 1]);

angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));

meanIntBins = linspace(0,1,21);
meanIntCents = meanIntBins+0.5*mean(diff(meanIntBins));
meanIntCents(end) = [];

for condID = 1:2
    for flyID = 1:cond{condID}.numFlies
        trialName = 'All_30C';
        meanIntAllDark = [];
        meanIntAllCL = [];
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))
            
            % Get the behavioral data
            tAll = cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OffsetRotMatch(:,1);
            tAll = tAll-tAll(1);
            heading = pi/180*cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.PosRotMatch;
            
            % Find the dark periods
            DarkPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Direction == 0);

            % Find the closed loop period
            CLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.Closed > 0);
            CLPer = setdiff(CLPer,DarkPer);
            CLPer(end) = [];

            % Find the rotational open loop period
            OLPer = find(cond{condID}.allFlyData{flyID}.(trialName){trialID}.positionDatMatch.OLGain > 0);
            
            % Get the activity
            EBAct = cond{condID}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax-1;
            
            % Calculate the PVA
            num_ROIs = size(EBAct,1);
            clear meanAngRaw;
            clear meanIntRaw;
            for ts = 1:length(EBAct)
                meanAngRaw(ts) = circ_mean(angs', squeeze(EBAct(:,ts)));
                meanIntRaw(ts) = circ_r(angs', squeeze(EBAct(:,ts)));
            end
            
            PVAPlt = meanAngRaw;
            for ts = 1:length(PVAPlt)-1
                if abs(PVAPlt(ts+1) - PVAPlt(ts)) > pi
                    PVAPlt(ts) = NaN;
                end
            end

            % Plot the PVA comparison in the dark
            figure(PVADark);
            subplot(5,3,trialID+(flyID-1)*3+(condID-1)*12)
            hold on;
            plot(tAll(DarkPer),PVAPlt(DarkPer),'Color',[0 0.75 0],'LineWidth',2);
            plot(tAll(DarkPer),heading(DarkPer),'k');
            xlim([tAll(DarkPer(1)) tAll(DarkPer(end))]);
            ylim([-pi pi]);
            if condID == 2
                xlabel('time (s)');
            end
            if trialID == 1
                ylabel('position (rad)');
            end
            if trialID == 1
                title(strcat(cond{condID}.name,' - fly #',num2str(flyID)))
                if flyID == 1 & condID == 1
                    text(tAll(DarkPer(1))-10,1.25*pi,'dark','FontSize',14);
                    legend({'PVA','heading'});
                    legend('boxoff');
                end
            end
            
            % Plot the PVA comparison with a stripe
            figure(PVACL);
            subplot(5,3,trialID+(flyID-1)*3+(condID-1)*12)
            hold on;
            plot(tAll(CLPer),PVAPlt(CLPer),'Color',[0 0.75 0],'LineWidth',2);
            plot(tAll(CLPer),heading(CLPer),'b');
            xlim([tAll(CLPer(1)) tAll(CLPer(end))]);
            ylim([-pi pi]);
            if condID == 2
                xlabel('time (s)');
            end
            if trialID == 1
                ylabel('position (rad)');
            end
            if trialID == 1
                title(strcat(cond{condID}.name,' - fly #',num2str(flyID)))
                if flyID == 1 & condID == 1
                    text(tAll(CLPer(1))-10,1.25*pi,'closed loop stripe','FontSize',14);
                    legend({'PVA','heading'});
                    legend('boxoff');
                end
            end
        end
    end
end