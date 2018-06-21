%% Clear out old data and load the new data
clear;
clc;

cond = FlyDatLoad(2);

%% Plot an example fly
condEx = 1;
flyEx = 3;
trialEx = 4;

tAll = cond{condEx}.allFlyData{flyEx}.All{trialEx}.positionDatMatch.OffsetRotMatch(:,1);

DarkPer = find(cond{condEx}.allFlyData{flyEx}.All{trialEx}.positionDatMatch.Direction == 0);

CLPer = find(cond{condEx}.allFlyData{flyEx}.All{trialEx}.positionDatMatch.Closed > 0);
CLPer = setdiff(CLPer,DarkPer);
CLPer(end) = [];    

tAll = tAll-tAll(DarkPer(1));
heading = cond{condEx}.allFlyData{flyEx}.All{trialEx}.positionDatMatch.OffsetRotMatch(:,2);
RAct = cond{condEx}.allFlyData{flyEx}.All{trialEx}.RROIaveMax;
GAct = cond{condEx}.allFlyData{flyEx}.All{trialEx}.GROIaveMax;
angs = linspace(-pi,pi,size(RAct,1)+1);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));

headingPlt = heading;
for tPt = 2:length(heading)
    if abs(headingPlt(tPt) - headingPlt(tPt-1)) > pi
        headingPlt(tPt-1) = NaN;
    end
end

% Plot the profiles
actEx = figure('units','normalized','outerposition',[0 0 1 1]);

% Plot the activity
subplot(4,1,1);
hold on;
plot(tAll(DarkPer),heading(DarkPer),'k');
ylabel('heading (rad)');
ylim([-pi pi]);
xlim([tAll(DarkPer(1)) tAll(DarkPer(end))]);
set(gca,'YTick',[-pi:pi/2:pi],'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});

subplot(4,1,2);
RIm = zeros([size(RAct) 3]);
RIm(:,:,1) = flipud((RAct-min(min(RAct)))./...
    (max(max(RAct))-min(min(RAct))));
RIm(:,:,3) = flipud((RAct-min(min(RAct)))./...
    (max(max(RAct))-min(min(RAct))));
image(tAll(DarkPer),angs,RIm(:,DarkPer,:));
title('G-E activity');
ylabel('EB position (rad)');
xlim([tAll(DarkPer(1)) tAll(DarkPer(end))]);
set(gca,'YTick',[-pi:pi/2:pi],'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});

subplot(4,1,3);
GIm = zeros([size(GAct) 3]);
GIm(:,:,2) = flipud((GAct-min(min(GAct)))./...
    (max(max(GAct))-min(min(GAct))));
image(tAll(DarkPer),angs,GIm(:,DarkPer,:));
title('E-PG activity');
ylabel('EB position (rad)');
xlim([tAll(DarkPer(1)) tAll(DarkPer(end))]);
set(gca,'YTick',[-pi:pi/2:pi],'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});


subplot(4,1,4);
overlayIm = RIm(:,DarkPer,:)+GIm(:,DarkPer,:);
image(tAll(DarkPer),angs,overlayIm);
xlim([tAll(DarkPer(1)) tAll(DarkPer(end))]);
title('combined activity');
ylabel('EB position (rad)');
xlabel('time (sec)');
colormap(brewermap(64, 'Blues'));
set(gca,'YTick',[-pi:pi/2:pi],'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});

set(actEx,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(actEx,'C:\Users\turnerevansd\Documents\RawAnalysis\G-E\TwoColorExample','-dpdf'); 

%% Sort the activity by the visual conditions

allAct = {};

for condID = 1:length(cond)
   allAct{condID}.name = cond{condID}.name;
   
   % Step through the flies
   for flyID = 1:cond{condID}.numFlies
       
      % Sort the data across trials
      for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)
         
         % Find the dark period
         darkPer = find(cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.Direction == 0);
         
         % Find the open loop period
         OLPer = find(cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.OLGain > 0);
         
         % Find the clockwise rotations
         CWPer = find(cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.Direction > 0);
         CWPer = intersect(OLPer,CWPer);
         % Find the counter-clockwise rotations
         CCWPer = find(cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.Direction < 0);
         CCWPer = intersect(OLPer,CCWPer);
         
         % Find the closed loop period
         CLPer = find(cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.Direction > 0);
         CLPer = setdiff(CLPer,CWPer);
         
         % Extract the behavioral parameters
         vR = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.vRot;
         vF = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.vF;
         tPts = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.OffsetRotMatch(:,1);
         stripePos = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.OffsetRotMatch(:,2);
         heading = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.PosRotMatch;
         
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
                 allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.RAct = ...
                     cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax(:,perNow);
                 allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.GAct = ...
                     cond{condID}.allFlyData{flyID}.All{trialID}.GROIaveMax(:,perNow);
         end
      end
   end
end

%% Look at the correlation between the PVA and bump amplitude above a given threshold

% Set a PVA amplitude threshold
PVAAmpThresh = 0.25;

% Set animal activity thresholds
minvF = 0.1;
minvR = pi/10;

% Calculate the angles for the ROIs
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + 0.5*mean(diff(angs));

% Make the figure;
PVACorr = figure('units','normalized','outerposition',[0 0 1 1]);

% Step through the conditions
for condID = 1:length(cond)
    
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        % Step through the periods
        for periodID = 1
        
            % Step through the trials
            for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

                % Get the behavioral data
                heading = pi/180*allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.heading;
                
                % Get the activity
                REBAct = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.RAct-1;
                GEBAct = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.GAct-1;
                
                % Calculate the PVA
                num_ROIs = size(REBAct,1);
                clear RPVA GPVA;
                clear RPVAAmp GPVAAmp;
                for ts = 1:length(REBAct)
                    RPVA(ts) = circ_mean(angs', squeeze(REBAct(:,ts)));
                    RPVAAmp(ts) = circ_r(angs', squeeze(REBAct(:,ts)));
                    GPVA(ts) = circ_mean(angs', squeeze(GEBAct(:,ts)));
                    GPVAAmp(ts) = circ_r(angs', squeeze(GEBAct(:,ts)));
                end

                % Keep the PVA constant if its amplitude falls below a certain
                % threshold
                RPVAConst = RPVA;
                RZeros = [];
                for ts = 2:length(RPVAConst)
                    if RPVAAmp(ts) < PVAAmpThresh
                        RPVAConst(ts) = RPVAConst(ts-1);
                        RZeros = vertcat(RZeros,ts);
                    end
                end
                GPVAConst = GPVA;
                GZeros = [];
                for ts = 2:length(GPVAConst)
                    if GPVAAmp(ts) < PVAAmpThresh
                        GPVAConst(ts) = GPVAConst(ts-1);
                        GZeros = vertcat(GZeros,ts);
                    end
                end
                
                allZeros = union(RZeros,GZeros);
                
%                 RPVAUnwrap = UnWrap(RPVAConst,2,0);
%                 RPVAUnwrap(allZeros) = [];
%                 GPVAUnwrap = UnWrap(GPVAConst,2,0);
%                 GPVAUnwrap(allZeros) = [];
                
                RPVAConst(allZeros) = [];
                GPVAConst(allZeros) = [];       
                
                CCs = corrcoef(RPVAConst,GPVAConst);
                subplot(2,3,1)
                hold on;
                scatter(5*(condID-1)+trialID,CCs(2,1),15,[0.75*(2-condID) 0.75*(condID-1) 0.75*(2-condID)],'filled');
                alpha(0.8);
                ylim([0 1]);
                xlim([0 11]);
                set(gca,'XTick',[1:10],'XTickLabel',{'1','2','3','4','5','6','7','8','9','10'});
                xlabel('fly #');
                ylabel('corr. coef.');
                title('PVA position');
                
                CCs = corrcoef(max(REBAct),max(GEBAct));
                subplot(2,3,2)
                hold on;
                scatter(5*(condID-1)+trialID,CCs(2,1),15,[0.75*(2-condID) 0.75*(condID-1) 0.75*(2-condID)],'filled');
                alpha(0.8);
                ylim([0 1]);
                xlim([0 11]);
                set(gca,'XTick',[1:10],'XTickLabel',{'1','2','3','4','5','6','7','8','9','10'});
                xlabel('fly #');
                ylabel('corr. coef.');
                title('bump amplitude');

            end
        end
    end
end
subplot(2,3,1);
scatter(4.5,0.1,20,[0.75 0 0.75],'filled');
alpha(0.8);
text(4.75,0.1,'= \color{magenta}G-E\color{black}, \color{green}E-PG');

subplot(2,3,1);
scatter(4.5,0.05,20,[0 0.75 0],'filled');
alpha(0.8);
text(4.75,0.05,'= \color{green}G-E\color{black}, \color{magenta}E-PG');

set(PVACorr,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(PVACorr,'C:\Users\turnerevansd\Documents\RawAnalysis\G-E\CorrCoeffs','-dpdf'); 