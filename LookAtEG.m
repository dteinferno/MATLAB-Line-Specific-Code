%% Clear out old data and load the new data
clear;
clc;

cond = FlyDatLoad(1);

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
             allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.act = ...
                 cond{condID}.allFlyData{flyID}.All{trialID}.ROIaveMax(:,perNow);
         end
      end
   end
end

%% Plot example trials
% Specify the EB example
EBFly = 4;
EBTrial = 5;

angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs+0.5*mean(diff(angs));

flyID = EBFly;
trialID = EBTrial;
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs+0.5*mean(diff(angs));
yLimits = [-pi pi];

for periodID = 1:3
    % Make a figure;
    actEx = figure('units','normalized','outerposition',[0 0 1 1]);

    if periodID < 3
        tPts = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.tPts;
        tPts = tPts - tPts(1);
        stripePos = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.stripePos;
        heading = pi/180*allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.heading;
    else
        tPts1 = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.tPts;
        tPts2 = allAct{1}.fly{flyID}.period{periodID+1}.trial{trialID}.tPts;
        tPts2 = tPts2 - tPts1(1);
        tPts1 = tPts1 - tPts1(1);
        tPts = vertcat(tPts1,tPts2);
        stripePos1 = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.stripePos;
        stripePos2 = allAct{1}.fly{flyID}.period{periodID+1}.trial{trialID}.stripePos;
        heading1 = pi/180*allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.heading;
        heading2 = pi/180*allAct{1}.fly{flyID}.period{periodID+1}.trial{trialID}.heading;
    end


    stripePos = -stripePos;
    heading = -heading;
    if periodID > 2
        stripePos1 = -stripePos1;
        stripePos2 = -stripePos2;
        heading1 = -heading1;
        heading2 = -heading2;
    end 

    subplot(4,1,1)
    hold on;
    if periodID < 3
        plot(tPts,stripePos,...
            'Color',allAct{1}.fly{flyID}.period{periodID}.color);
    else
        plot(tPts1,stripePos1,...
            'Color',allAct{1}.fly{flyID}.period{periodID}.color);
        plot(tPts2,stripePos2,...
            'Color',allAct{1}.fly{flyID}.period{periodID+1}.color);
        plot(tPts1,heading1,...
            'Color','k');
        plot(tPts2,heading2,...
            'Color','k');
    end
    ylim([-pi pi]);
    xlim([tPts(1) tPts(end)]);
    ylabel('position (rad)');
    title(allAct{1}.fly{flyID}.period{periodID}.type);
    colorbar;

    subplot(4,1,2)
    Act = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.act-1;
    if periodID ==3 
        Act2 = allAct{1}.fly{flyID}.period{periodID+1}.trial{trialID}.act-1;
        Act = horzcat(Act,Act2);
    end

    imagesc(tPts,angs,Act);
    colorbar;
    caxis([0 1.5]);
    ylim(yLimits);
    xlim([tPts(1) tPts(end)]);
    ylabel('EB position');
    xlabel('time (s)');
    
    colormap(brewermap(64, 'Blues'));
    
     set(actEx,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
        print(actEx,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\',...
            allAct{1}.name,'\',allAct{1}.name,'_',...
            allAct{1}.fly{flyID}.period{periodID}.type,'_Ex'),'-dpdf');
end

%% Run a regression with the peak maximum against the forward and rotational velocities - for each trial
numPts = 41;

% Look at the EB and PB
for condID = 1:length(cond)
   reg = figure('units','normalized','outerposition',[0 0 1 1]);
   
   numFlies = length(allAct{condID}.fly);
   % Step through the flies
   for flyID = 1:cond{condID}.numFlies
       
      % Sort the data across trials
      for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)
          
          numPeriods = length(allAct{condID}.fly{flyID}.period);
          
          % Sort the data by period
          for periodID = 1:4
              
                  % Find the time step
                  tStep = mean(diff(allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.tPts));
                  
                  % Get the variables of choice
                  vR = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vR;
                  vF = allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vF;
                  
                  % Get the maximum activity
                  maxAct = max(allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.act);

                  autoC = zeros(numPts,1);
                  for lag = 1:numPts
                      lagCCs = corrcoef(...
                          abs(vR(ceil(numPts/2):end-floor(numPts/2)-1)),...
                          maxAct(numPts-lag+1:end-lag));
                      autoC(lag) = lagCCs(2,1);
                  end
                  allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vRAutoC = autoC;

                  subplot(2*numFlies,numPeriods,2*numPeriods*(flyID-1)+periodID);
                  hold on;
                  plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoC,'Color',...
                      allAct{condID}.fly{flyID}.period{periodID}.color);
                  line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                  line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                  ylim([-0.5 1]);
                  xlim([-tStep*numPts/2 tStep*numPts/2]);
                  if periodID == 1
                      ylabel('|vR| autocorrelation');
                      text(-2*tStep*numPts/2,1.1,strcat('fly #',num2str(flyID)));
                      if flyID == 1
                          text(-tStep*numPts/2,1.25,allAct{condID}.name);
                      end
                  end
                  if flyID == 1
                      title(strcat(allAct{condID}.fly{flyID}.period{periodID}.type,'-red'),...
                          'Color',allAct{condID}.fly{flyID}.period{periodID}.color);
                  end

                  autoC = zeros(numPts,1);
                  for lag = 1:numPts
                      lagCCs = corrcoef(...
                          vF(ceil(numPts/2):end-floor(numPts/2)-1),...
                          maxAct(numPts-lag+1:end-lag));
                      autoC(lag) = lagCCs(2,1);
                  end
                  allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vFAutoC = autoC;

                  subplot(2*numFlies,numPeriods,2*numPeriods*(flyID-1)+periodID+numPeriods);
                  hold on;
                  plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoC,'Color',...
                      allAct{condID}.fly{flyID}.period{periodID}.color);
                  line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                  line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                  ylim([-0.5 1]);
                  xlim([-tStep*numPts/2 tStep*numPts/2]);
                  if periodID == 1
                      ylabel('vF autocorrelation');
                  end
                  if flyID == cond{condID}.numFlies
                      xlabel('time delay (s)');
                  end
              end
          end
      end
end
          
%% Run a regression with the peak maximum against the forward and rotational velocities - averages

% Look at the EB and PB
for condID = 1:length(cond)
   reg = figure('units','normalized','outerposition',[0 0 1 1]);
   
   numFlies = length(allAct{condID}.fly);
   % Step through the flies
   for flyID = 1:cond{condID}.numFlies
      % Sort the data by period
      for periodID = 1:4

              vRAll = [];
              vFAll = [];
                  
              % Sort the data across trials
              for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

                  % Group across trials
                  vRAll = horzcat(vRAll,allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vRAutoC);
                  vFAll = horzcat(vFAll,allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vFAutoC);
                  
              end
               
                subplot(2,numFlies,flyID);
                hold on;
                plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRAll,2),'Color',...
                allAct{condID}.fly{flyID}.period{periodID}.color);
            if flyID == 1
                    ylabel('|vR| autocorrelation');
            end
                if flyID == 1
                    text(-tStep*numPts/2,1.25,allAct{condID}.name);
                    if periodID == 4
                       legend({'dark','CL','CW','CCW'});
                    end
                end
                    title(strcat('fly #',num2str(flyID),'-red'));
                if periodID == 4
                    line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                    line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                    ylim([-0.25 0.5]);
                    xlim([-tStep*numPts/2 tStep*numPts/2]);
                end

                subplot(2,numFlies,flyID+numFlies);
                hold on;
                plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFAll,2),'Color',...
                  allAct{condID}.fly{flyID}.period{periodID}.color);
                line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                ylim([-0.25 0.5]);
                xlim([-tStep*numPts/2 tStep*numPts/2]);
                if  flyID == 1
                    ylabel('vF autocorrelation');
                end
                xlabel('time delay (s)');
      end
  end
end   

%% In the EB, look at the mean peak max vs. vRot - binned

% Specify rotational velocity bin edges and centers
vREdges = linspace(0.001,pi/2,9);
vRCents = vREdges(1:end-1)+0.5*mean(diff(vREdges));
vREdges = horzcat(0,vREdges);
vRCents = horzcat(0.0005,vRCents);


% Specify Savtisky-Golay filter parameters for smoothing the velocity
sgolayOrder = 3;
sgolayFrames = 11;

% Look at the EB condition specifically
for condID = 1
    actMeanAll = zeros(4,length(vRCents));
    actSDAll = zeros(4,length(vRCents));
            
    actBin = figure('units','normalized','outerposition',[0 0 1 1]);
    numFlies = length(allAct{condID}.fly);
    for flyID = 1:numFlies
        for periodID = 1:4
            % Create arrays to hold max velocities and bump amplitudes
            vRAll = [];
            maxActAll = [];
            for trialID = 1:length(allAct{condID}.fly{flyID}.period{periodID}.trial)
                vR = abs(sgolayfilt(allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vR,...
                    sgolayOrder,sgolayFrames));
                vRAll = vertcat(vRAll,vR);

                maxAct = max(allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.act)-1;
                maxActAll = vertcat(maxActAll,maxAct'); 
            end
            binIDs = discretize(vRAll,vREdges);
            for binNow = 1:length(vREdges)-1
                meanAct = mean(maxActAll(find(binIDs == binNow)));
                actMeanAll(periodID,binNow) = meanAct;
                SDAct = std(maxActAll(find(binIDs == binNow)));
                actSDAll(periodID,binNow) = SDAct;
            end
        end

        subplot(3,numFlies,flyID);
        hold on;
        title(strcat('fly #',num2str(flyID)));
        for periodID = 1:4
            meanNow = squeeze(actMeanAll(periodID,:));
            SDNow = squeeze(actSDAll(periodID,:));
            plot(vRCents,meanNow,...
                'color',allAct{condID}.fly{flyID}.period{periodID}.color);
            p = patch('XData',[vRCents fliplr(vRCents)],...
                'YData',[meanNow+SDNow fliplr(meanNow-SDNow)],...
                'FaceColor',allAct{condID}.fly{flyID}.period{periodID}.color,...
                'FaceAlpha',0.1,...
                'EdgeColor','none');
            uistack(p,'bottom');

        end
        xlim([vREdges(1) vREdges(end)]);
        ylim([0 1.25]);
        line([0 0],[0 1.25],'Color','k','LineStyle','--');
                
        if flyID == 1
            subplot(3,numFlies,flyID)
            ylabel('max DF/F');
        end
    end
end

%% In the PB, look at the mean peak max vs. vF - binned

% Specify rotational velocity bin edges and centers
vFEdges = linspace(0.001,0.5,9);
vFCents = vFEdges(1:end-1)+0.5*mean(diff(vFEdges));
vFEdges = horzcat(0,vFEdges);
vFCents = horzcat(0,vFCents);

% Specify Savtisky-Golay filter parameters for smoothing the velocity
sgolayOrder = 3;
sgolayFrames = 11;

% Look at the EB condition specifically
for condID = 1
    actMeanAll = zeros(4,length(vRCents));
    actSDAll = zeros(4,length(vRCents));
            
    actBin = figure('units','normalized','outerposition',[0 0 1 1]);
    numFlies = length(allAct{condID}.fly);
    for flyID = 1:numFlies
        for periodID = 1:4
            % Create arrays to hold max velocities and bump amplitudes
            vFAll = [];
            maxActAll = [];
            for trialID = 1:length(allAct{condID}.fly{flyID}.period{periodID}.trial)
                vF = sgolayfilt(allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.vF,...
                    sgolayOrder,sgolayFrames);
                vFAll = vertcat(vFAll,vF);

                maxAct = max(allAct{condID}.fly{flyID}.period{periodID}.trial{trialID}.act)-1;
                maxActAll = vertcat(maxActAll,maxAct'); 
            end
            binIDs = discretize(vFAll,vFEdges);
            for binNow = 1:length(vFEdges)-1
                meanAct = mean(maxActAll(find(binIDs == binNow)));
                actMeanAll(periodID,binNow) = meanAct;
                SDAct = std(maxActAll(find(binIDs == binNow)));
                actSDAll(periodID,binNow) = SDAct;
            end
        end

        subplot(3,numFlies,flyID);
        hold on;
        title(strcat('fly #',num2str(flyID)));
        for periodID = 1:4
            meanNow = squeeze(actMeanAll(periodID,:));
            SDNow = squeeze(actSDAll(periodID,:));
            plot(vFCents,meanNow,...
                'color',allAct{condID}.fly{flyID}.period{periodID}.color);
            p = patch('XData',[vFCents fliplr(vFCents)],...
                'YData',[meanNow+SDNow fliplr(meanNow-SDNow)],...
                'FaceColor',allAct{condID}.fly{flyID}.period{periodID}.color,...
                'FaceAlpha',0.1,...
                'EdgeColor','none');
            uistack(p,'bottom');

        end
        xlim([vFEdges(1) vFEdges(end)]);
        ylim([0 1.25]);
        line([0 0],[0 1.25],'Color','k','LineStyle','--');
                
        if flyID == 1
            subplot(3,numFlies,flyID)
            ylabel('max DF/F');
        end
    end
end

%% Look at the bump persistance


% Set the Savitsky Golay filter parameters
sgolayOrder = 3;
sgolayFrames = 11;

% Set the velocity thresholds
vRThresh = 0.01*pi;
vFThresh = 0.01;

% Set a threshold for the length of a standing bout
minBoutL = 5;

% Step through the flies
for flyID = 1:cond{1}.numFlies

    % Sort the data by period
    for periodID = 1:4
        
        numPeriods = length(allAct{1}.fly{flyID}.period);
        maxBoutL = 0; % Max bout length
        numBouts = 0; % Number of bouts
            
        % Sort the data across trials
        for trialID = 1:length(cond{1}.allFlyData{flyID}.All)

            % Find the time step
            tPts = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.tPts;
            tStep = mean(diff(tPts));

            % Find when the animal is standing
            vR = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.vR;
            vR = sgolayfilt(vR,sgolayOrder,sgolayFrames);
            vF = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.vF;
            vF = sgolayfilt(vF,sgolayOrder,sgolayFrames);
            standingBouts = intersect(find(abs(vR)<=vRThresh),find(vF<=vFThresh));

            % Group the standing bouts
            boutStart = vertcat(1,find(diff(standingBouts)>1)+1);
            boutStop = vertcat(find(diff(standingBouts)>1),length(standingBouts));

            % Remove shorter bouts
            shortBouts = [];
            for bout = 1:length(boutStart)
             if (boutStop(bout)-boutStart(bout)) <= minBoutL 
                 shortBouts = vertcat(shortBouts,bout);
             end
            end
            boutStart(shortBouts) = [];
            boutStop(shortBouts) = [];
            
            maxBoutL = max(maxBoutL,max(boutStop-boutStart));
            numBouts = numBouts+length(boutStart);

            % Add the standing bout info to the object
            allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.boutStart = standingBouts(boutStart);
            allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.boutStop = standingBouts(boutStop);

        end
        
        allAct{1}.fly{flyID}.period{periodID}.maxBoutL = maxBoutL;
        allAct{1}.fly{flyID}.period{periodID}.numBouts = numBouts;
        
    end
end

% Group the standing bouts across trials

% Make the figure
bumpPer = figure('units','normalized','outerposition',[0 0 1 1]);

% Step through the flies
for flyID = 1:cond{1}.numFlies
    
    % Sort the data by period
    for periodID = 1:4
    
        numPeriods = length(allAct{1}.fly{flyID}.period);
        maxBoutL = allAct{1}.fly{flyID}.period{periodID}.maxBoutL;
        numBouts = allAct{1}.fly{flyID}.period{periodID}.numBouts;
        
        standAct = zeros(numBouts,16,maxBoutL);
        standActMean = zeros(16,maxBoutL);
        standBout = 1;

        % Sort the data across trials
        for trialID = 1:length(cond{1}.allFlyData{flyID}.All)

            % Place all of the standing bout activity into the array
            actNow = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.act-1;
            boutStart = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.boutStart;
            boutStop = allAct{1}.fly{flyID}.period{periodID}.trial{trialID}.boutStop;
            
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
        
        subplot(numPeriods,cond{1}.numFlies,flyID+cond{1}.numFlies*(periodID-1))
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
        title(allAct{1}.fly{flyID}.period{periodID}.type);
    end
end

set(bumpPer,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpPer,'C:\Users\turnerevansd\Documents\RawAnalysis\G-E\BumpPersistance','-dpdf');