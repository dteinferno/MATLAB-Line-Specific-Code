%% Clear out old data and load the new data
clear;
clc;

cond = FlyDatLoad(2);

%% Specify the PB glomeruli of interest
RLPB = [2:9];
RRPB = [10:17];
GLPB = [1:8];
GRPB = [11:18];

%% Sort the activity by the visual conditions

allAct = {};
% Look at the EB and PB
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
             
             for colorID = 1:2
                 if periodID == 1
                     perNow = darkPer;
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type = 'dark';
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color = [0 0 0];
                 elseif periodID == 2
                     perNow = CLPer;
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type = 'CL';
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color = [0 0 1];
                 elseif periodID == 3
                     perNow = CWPer;
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type = 'CW';
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color = [0.5 0 1];
                 elseif periodID == 4
                     perNow = CCWPer;
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type = 'CCW';
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color = [0 0.5 1];
                 end
                 allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vR = vR(perNow);
                 allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vF = vF(perNow);
                 
                 % Find the time points, stripe position, and heading
                 allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.tPts = tPts(perNow);
                 allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.stripePos = stripePos(perNow);
                 allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.heading = heading(perNow);
             end
             % Find the activity
             if contains(cond{condID}.name,'EB')
                 allAct{condID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.act = ...
                     cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax(:,perNow);
                 allAct{condID}.fly{flyID}.color{2}.period{periodID}.trial{trialID}.act = ...
                     cond{condID}.allFlyData{flyID}.All{trialID}.GROIaveMax(:,perNow);
             elseif contains(cond{condID}.name,'PB')
                 for RL = 1:2
                     if RL == 1
                         RROI = RLPB;
                         GROI = GLPB;
                     else
                         RROI = RRPB;
                         GROI = GRPB;
                     end
                     allAct{condID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.PBSide{RL}.act =...
                         cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax(RROI,perNow);
                     allAct{condID}.fly{flyID}.color{2}.period{periodID}.trial{trialID}.PBSide{RL}.act =...
                         cond{condID}.allFlyData{flyID}.All{trialID}.GROIaveMax(GROI,perNow);
                 end
             end
         end
      end
   end
end

%% Plot example trials
% Specify the EB example
EBFly = 2;
EBTrial = 2;

% Specify the PB example
PBFly = 3;
PBTrial = 3;

angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs+0.5*mean(diff(angs));

% Plot the dark periods
for regionID = 1:2
    % Specify EB and PB specifics
    if regionID == 1
        flyID = EBFly;
        trialID = EBTrial;
        angs = linspace(-pi,pi,17);
        angs(end) = [];
        angs = angs+0.5*mean(diff(angs));
        yLimits = [-pi pi];
    else
        flyID = PBFly;
        trialID = PBTrial;
        angs = [1:18];
        yLimits = [1 18];
    end
    for periodID = 1:3
        % Make a figure;
        if periodID < 4
            actEx = figure('units','normalized','outerposition',[0 0 1 1]);
        end
    
        if periodID < 3
            tPts = allAct{regionID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.tPts;
            tPts = tPts - tPts(1);
            stripePos = allAct{regionID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.stripePos;
            heading = pi/180*allAct{regionID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.heading;
        else
            tPts1 = allAct{regionID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.tPts;
            tPts2 = allAct{regionID}.fly{flyID}.color{1}.period{periodID+1}.trial{trialID}.tPts;
            tPts2 = tPts2 - tPts1(1);
            tPts1 = tPts1 - tPts1(1);
            tPts = vertcat(tPts1,tPts2);
            stripePos1 = allAct{regionID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.stripePos;
            stripePos2 = allAct{regionID}.fly{flyID}.color{1}.period{periodID+1}.trial{trialID}.stripePos;
            heading1 = pi/180*allAct{regionID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.heading;
            heading2 = pi/180*allAct{regionID}.fly{flyID}.color{1}.period{periodID+1}.trial{trialID}.heading;
        end
        
        
        if regionID == 1
            stripePos = -stripePos;
            heading = -heading;
            if periodID > 2
                stripePos1 = -stripePos1;
                stripePos2 = -stripePos2;
                heading1 = -heading1;
                heading2 = -heading2;
            end 
        end
        
        subplot(4,1,1)
        hold on;
        if periodID < 3
            plot(tPts,stripePos,...
                'Color',allAct{regionID}.fly{flyID}.color{1}.period{periodID}.color);
        else
            plot(tPts1,stripePos1,...
                'Color',allAct{regionID}.fly{flyID}.color{1}.period{periodID}.color);
            plot(tPts2,stripePos2,...
                'Color',allAct{regionID}.fly{flyID}.color{1}.period{periodID+1}.color);
            plot(tPts1,heading1,...
                'Color','k');
            plot(tPts2,heading2,...
                'Color','k');
        end
        ylim([-pi pi]);
        xlim([tPts(1) tPts(end)]);
        ylabel('position (rad)');
        title(allAct{regionID}.fly{flyID}.color{1}.period{periodID}.type);
        colorbar;

        subplot(4,1,2)
        if regionID == 1
            RAct = allAct{regionID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.act-1;
        else
            RActLPB = allAct{regionID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.PBSide{1}.act-1;
            RActRPB = allAct{regionID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.PBSide{2}.act-1;
            RAct = zeros([20, length(RActLPB)]);
            RAct(RLPB,:) = RActLPB;
            RAct(RRPB+2,:) = RActRPB;
        end
        if periodID ==3 
           if regionID == 1
                RAct2 = allAct{regionID}.fly{flyID}.color{1}.period{periodID+1}.trial{trialID}.act-1;
                RAct = horzcat(RAct,RAct2);
            else
                RActLPB2 = allAct{regionID}.fly{flyID}.color{1}.period{periodID+1}.trial{trialID}.PBSide{1}.act-1;
                RActRPB2 = allAct{regionID}.fly{flyID}.color{1}.period{periodID+1}.trial{trialID}.PBSide{2}.act-1;
                RAct2 = zeros([20, length(RActLPB2)]);
                RAct2(RLPB,:) = RActLPB2;
                RAct2(RRPB+2,:) = RActRPB2;
                RAct = horzcat(RAct,RAct2);
            end 
        end
            
        imagesc(tPts,angs,RAct);
        colorbar;
        caxis([0 1]);
        ylim(yLimits);
        xlim([tPts(1) tPts(end)]);
        if regionID == 1
            ylabel('EB position');
        else
            ylabel('R PB        L PB');
        end
        title('P-EG act.');

        subplot(4,1,3)
        if regionID == 1
            GAct = allAct{regionID}.fly{flyID}.color{2}.period{periodID}.trial{trialID}.act-1;
        else
            GActLPB = allAct{regionID}.fly{flyID}.color{2}.period{periodID}.trial{trialID}.PBSide{1}.act-1;
            GActRPB = allAct{regionID}.fly{flyID}.color{2}.period{periodID}.trial{trialID}.PBSide{2}.act-1;
            GAct = zeros([20, length(GActLPB)]);
            GAct(GLPB,:) = GActLPB;
            GAct(GRPB+2,:) = GActRPB;
        end
        if periodID ==3
           if regionID == 1
                GAct2 = allAct{regionID}.fly{flyID}.color{2}.period{periodID+1}.trial{trialID}.act-1;
                GAct = horzcat(GAct,GAct2);
            else
                GActLPB2 = allAct{regionID}.fly{flyID}.color{2}.period{periodID+1}.trial{trialID}.PBSide{1}.act-1;
                GActRPB2 = allAct{regionID}.fly{flyID}.color{2}.period{periodID+1}.trial{trialID}.PBSide{2}.act-1;
                GAct2 = zeros([20, length(GActLPB2)]);
                GAct2(GLPB,:) = GActLPB2;
                GAct2(GRPB+2,:) = GActRPB2;
                GAct = horzcat(GAct,GAct2);
            end 
        end
        imagesc(tPts,angs,GAct);
        colorbar;
        caxis([0 1]);
        ylim(yLimits);
        xlim([tPts(1) tPts(end)]);
        if regionID == 1
            ylabel('EB position');
        else
            ylabel('R PB        L PB');
        end
        title('P-EN2 act.');

        colormap(brewermap(64, 'Blues'));

        subplot(4,1,4)
        overlayIm = zeros([size(RAct) 3]);
        overlayIm(:,:,1) = RAct;
        overlayIm(:,:,3) = RAct;
        overlayIm(:,:,2) = GAct;
        image(tPts,angs,overlayIm);
        colorbar;
        ylim(yLimits);
        xlim([tPts(1) tPts(end)]);
        if regionID == 1
            ylabel('EB position');
        else
            ylabel('R PB        L PB');
        end
        xlabel('time (s)');
        
        set(actEx,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
        print(actEx,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\',...
        allAct{regionID}.name,'\',allAct{regionID}.name,'_',...
        allAct{regionID}.fly{flyID}.color{1}.period{periodID}.type,'_Ex'),'-dpdf');
    end
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
          
          numPeriods = length(allAct{condID}.fly{flyID}.color{1}.period);
          
          % Sort the data by period
          for periodID = 1:4
              
              % Sort the data by color
              for colorID = 1:2
                  
                  % Find the time step
                  tStep = mean(diff(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.tPts));
                  
                  % Get the variables of choice
                  vR = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vR;
                  vF = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vF;
                  
                  % Get the maximum activity
                  if contains(allAct{condID}.name,'EB')
                      maxAct = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.act);
                      
                      autoC = zeros(numPts,1);
                      for lag = 1:numPts
                          lagCCs = corrcoef(...
                              abs(vR(ceil(numPts/2):end-floor(numPts/2)-1)),...
                              maxAct(numPts-lag+1:end-lag));
                          autoC(lag) = lagCCs(2,1);
                      end
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vRAutoC = autoC;
                      
                      subplot(2*numFlies,2*numPeriods,4*numPeriods*(flyID-1)+2*periodID-2+colorID);
                      hold on;
                      plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoC,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                      line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                      ylim([-0.5 1]);
                      xlim([-tStep*numPts/2 tStep*numPts/2]);
                      if periodID == 1 & colorID == 1
                          ylabel('|vR| autocorrelation');
                          text(-2*tStep*numPts/2,1.1,strcat('fly #',num2str(flyID)));
                          if flyID == 1
                              text(-tStep*numPts/2,1.25,allAct{condID}.name);
                          end
                      end
                      if flyID == 1
                          if colorID == 1
                              title(strcat(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type,'-red'),...
                                  'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                          else
                              title(strcat(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type,'-green'),...
                                  'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                          end
                      end
                      
                      autoC = zeros(numPts,1);
                      for lag = 1:numPts
                          lagCCs = corrcoef(...
                              vF(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxAct(numPts-lag+1:end-lag));
                          autoC(lag) = lagCCs(2,1);
                      end
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vFAutoC = autoC;
                      
                      subplot(2*numFlies,2*numPeriods,4*numPeriods*(flyID-1)+2*periodID-2+colorID+2*numPeriods);
                      hold on;
                      plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoC,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                      line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                      ylim([-0.5 1]);
                      xlim([-tStep*numPts/2 tStep*numPts/2]);
                      if periodID == 1 & colorID == 1
                          ylabel('vF autocorrelation');
                      end
                      if flyID == cond{condID}.numFlies
                          xlabel('time delay (s)');
                      end
                      
                      
                  elseif contains(allAct{condID}.name,'PB')
                      maxActL = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.act);
                      maxActR = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.act);
                      
                      vRPos = zeros(size(vR));
                      vRPos(find(vR>0)) = vR(find(vR > 0));
                      autoCL = zeros(numPts,1);
                      autoCR = zeros(numPts,1);
                      for lag = 1:numPts
                          
                          lagCCsL = corrcoef(...
                              vRPos(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActL(numPts-lag+1:end-lag));
                          autoCL(lag) = lagCCsL(2,1);
                          lagCCsR = corrcoef(...
                              vRPos(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActR(numPts-lag+1:end-lag));
                          autoCR(lag) = lagCCsR(2,1);
                      end
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRPosAutoC = autoCL;
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRPosAutoC = autoCR;
                      
                      subplot(3*numFlies,2*numPeriods,6*numPeriods*(flyID-1)+2*periodID-2+colorID);
                      hold on;
                      Lplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCL,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      Rplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCR,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                          'LineWidth',2,'LineStyle','--');
                      line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                      line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                      ylim([-0.5 1]);
                      xlim([-tStep*numPts/2 tStep*numPts/2]);
                      if periodID == 1 & colorID == 1
                          ylabel('vR > 0 autocorrelation');
                            text(-2*tStep*numPts/2,1.1,strcat('fly #',num2str(flyID)));
                            if flyID == 1
                                text(-tStep*numPts/2,1.25,allAct{condID}.name);
                                legend([Lplt Rplt],{'L PB','R PB'});
                            end
                      end
                      if flyID == 1
                          if colorID == 1
                              title(strcat(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type,'-red'),...
                                  'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                          else
                              title(strcat(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type,'-green'),...
                                  'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                          end
                      end
                      
                      vRNeg = zeros(size(vR));
                      vRNeg(find(vR < 0)) = -vR(find(vR < 0));
                      autoCL = zeros(numPts,1);
                      autoCR = zeros(numPts,1);
                      for lag = 1:numPts
                          
                          lagCCsL = corrcoef(...
                              vRNeg(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActL(numPts-lag+1:end-lag));
                          autoCL(lag) = lagCCsL(2,1);
                          lagCCsR = corrcoef(...
                              vRNeg(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActR(numPts-lag+1:end-lag));
                          autoCR(lag) = lagCCsR(2,1);
                      end
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRNegAutoC = autoCL;
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRNegAutoC = autoCR;
                      
                      subplot(3*numFlies,2*numPeriods,6*numPeriods*(flyID-1)+2*periodID-2+colorID+2*numPeriods);
                      hold on;
                      Lplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCL,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      Rplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCR,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                          'LineWidth',2,'LineStyle','--');
                      line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                      line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                      ylim([-0.5 1]);
                      xlim([-tStep*numPts/2 tStep*numPts/2]);
                      if periodID == 1 & colorID == 1
                          ylabel('vR < 0 autocorrelation');
                      end
                      
                      autoCL = zeros(numPts,1);
                      autoCR = zeros(numPts,1);
                      for lag = 1:numPts
                          lagCCsL = corrcoef(...
                              vF(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActL(numPts-lag+1:end-lag));
                          autoCL(lag) = lagCCsL(2,1);
                          lagCCsR = corrcoef(...
                              vF(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActR(numPts-lag+1:end-lag));
                          autoCR(lag) = lagCCsR(2,1);
                      end
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vFAutoC = autoCL;
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vFAutoC = autoCR;
                      
                      subplot(3*numFlies,2*numPeriods,6*numPeriods*(flyID-1)+2*periodID-2+colorID+4*numPeriods);
                      hold on;
                      Lplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCL,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      Rplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCR,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                          'LineWidth',2,'LineStyle','--');
                      line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                      line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                      ylim([-0.5 1]);
                      xlim([-tStep*numPts/2 tStep*numPts/2]);
                      if periodID == 1 & colorID == 1
                          ylabel('vF autocorrelation');
                      end
                      if flyID == cond{condID}.numFlies
                          xlabel('time delay (s)');
                      end
                  end
                  
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

          % Sort the data by color
          for colorID = 1:2
              if contains(allAct{condID}.name,'EB')
                  vRAll = [];
                  vFAll = [];
              elseif contains(allAct{condID}.name,'PB')
                  vRPosLAll = [];
                  vRNegLAll = [];
                  vFLAll = [];
                  vRPosRAll = [];
                  vRNegRAll = [];
                  vFRAll = [];
              end
              % Sort the data across trials
              for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

                  % Group across trials
                  if contains(allAct{condID}.name,'EB')
                      vRAll = horzcat(vRAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vRAutoC);
                      vFAll = horzcat(vFAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vFAutoC);
                  elseif contains(allAct{condID}.name,'PB') 
                      vRPosLAll = horzcat(vRPosLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRPosAutoC);
                      vRPosRAll = horzcat(vRPosLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRPosAutoC);

                      vRNegLAll = horzcat(vRNegLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRNegAutoC);
                      vRNegRAll = horzcat(vRNegLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRNegAutoC);

                      vFLAll = horzcat(vFLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vFAutoC);
                      vFRAll = horzcat(vFRAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vFAutoC);
                  end
              end
              if contains(allAct{condID}.name,'EB')
                    subplot(2,2*numFlies,2*(flyID-1)+colorID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRAll,2),'Color',...
                    allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    if colorID == 1
                        ylabel('|vR| autocorrelation');
                        if flyID == 1
                            text(-tStep*numPts/2,0.55,allAct{condID}.name);
                            if periodID == 4
                               legend({'dark','CL','CW','CCW'});
                            end
                        end
                        title(strcat('fly #',num2str(flyID),'-red'));
                    else
                        title(strcat('fly #',num2str(flyID),'-green'));
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 0.5]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                    end
                    
                    subplot(2,2*numFlies,2*(flyID-1)+colorID+2*numFlies);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                    line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                    ylim([-0.25 0.5]);
                    xlim([-tStep*numPts/2 tStep*numPts/2]);
                    if colorID == 1 & flyID == 1
                        ylabel('vF autocorrelation');
                    end
                    xlabel('time delay (s)');
              elseif contains(allAct{condID}.name,'PB')
                    subplot(3,2*numFlies,2*(flyID-1)+colorID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRPosRAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRPosLAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    if colorID == 1
                        ylabel('vR > 0 autocorrelation');
                        if flyID == 1
                            text(-tStep*numPts/2,0.6,allAct{condID}.name);
                        end
                        title(strcat('fly #',num2str(flyID),'-red'));
                    else
                        title(strcat('fly #',num2str(flyID),'-green'));
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 0.5]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                    end
                    
                    subplot(3,2*numFlies,2*(flyID-1)+colorID+2*numFlies);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRNegRAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRNegLAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    if colorID == 1
                        ylabel('vR < 0 autocorrelation');
                        title(strcat('fly #',num2str(flyID),'-red'));
                    else
                        title(strcat('fly #',num2str(flyID),'-green'));
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 0.5]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                    end

                    subplot(3,2*numFlies,2*(flyID-1)+colorID+4*numFlies);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFRAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFLAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                    line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                    ylim([-0.25 0.5]);
                    xlim([-tStep*numPts/2 tStep*numPts/2]);
                    if colorID == 1 & flyID == 1
                        ylabel('vF autocorrelation');
                    end
                    xlabel('time delay (s)');
              end
          end
      end
   end   
   set(reg,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    print(reg,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\',...
    allAct{condID}.name,'\',allAct{condID}.name,'_MeanAutoCC'),'-dpdf');
end

%% Plot the regressions nicely
% For each condition, pick and plot an example fly, and the correlation
% offset and delay for each condition

EBExFly = 5;
PBExFly = 2;

% Look at the EB and PB
for condID = 1:length(cond)
   reg = figure('units','normalized','outerposition',[0 0 1 1]);
   
   numFlies = length(allAct{condID}.fly);
   % Step through the flies
   for flyID = 1:cond{condID}.numFlies
      % Sort the data by period
      for periodID = 1:4

          % Sort the data by color
          for colorID = 1:2
              if contains(allAct{condID}.name,'EB')
                  vRAll = [];
                  vFAll = [];
              elseif contains(allAct{condID}.name,'PB')
                  vRPosLAll = [];
                  vRNegLAll = [];
                  vFLAll = [];
                  vRPosRAll = [];
                  vRNegRAll = [];
                  vFRAll = [];
              end
              % Sort the data across trials
              for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

                  % Group across trials
                  if contains(allAct{condID}.name,'EB')
                      vRAll = horzcat(vRAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vRAutoC);
                      vFAll = horzcat(vFAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vFAutoC);
                  elseif contains(allAct{condID}.name,'PB') 
                      vRPosLAll = horzcat(vRPosLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRPosAutoC);
                      vRPosRAll = horzcat(vRPosLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRPosAutoC);

                      vRNegLAll = horzcat(vRNegLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRNegAutoC);
                      vRNegRAll = horzcat(vRNegLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRNegAutoC);

                      vFLAll = horzcat(vFLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vFAutoC);
                      vFRAll = horzcat(vFRAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vFAutoC);
                  end
              end
              if contains(allAct{condID}.name,'EB')
                  if flyID == EBExFly
                      subplot(2,4,colorID);
                      hold on;
                      plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRAll,2),'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      ax = gca;
                      if colorID == 1
                          ylabel('|vR| autocorrelation');
                          text(-tStep*numPts/2-1,1.15,strcat(allAct{condID}.name,'-Summary'),'FontSize',14);
                          if periodID == 4
                             legend({'dark','CL','CW','CCW'});
                             legend('boxoff');
                          end
                          title('red channel');
                          ax.XColor = [0.5 0 0.5];
                          ax.YColor = [0.5 0 0.5];
                      else
                          title('green channel');
                          ax.XColor = [0 0.5 0];
                          ax.YColor = [0 0.5 0];
                      end
                      if periodID == 4
                          line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                          line([0 0],[-0.25 1], 'Color', 'k','LineStyle','--');
                          ylim([-0.25 1]);
                          xlim([-tStep*numPts/2 tStep*numPts/2]);
                      end
                                          
                      subplot(2,4,4+colorID);
                      hold on;
                      plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      if periodID == 4
                          line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                          line([0 0],[-0.25 1], 'Color', 'k','LineStyle','--');
                          ylim([-0.25 1]);
                          xlim([-tStep*numPts/2 tStep*numPts/2]);
                      end
                      ax = gca;
                      if colorID == 1
                          ylabel('vF autocorrelation');
                          ax.XColor = [0.5 0 0.5];
                          ax.YColor = [0.5 0 0.5];
                      else
                          ax.XColor = [0 0.5 0];
                          ax.YColor = [0 0.5 0];
                      end
                      xlabel('time delay (s)');
                  end
                           
                  for vVar = 1:2
                      if vVar ==1 
                          maxCC = max(mean(vRAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vRAll,2)==maxCC)-1)*tStep;
                      else
                          maxCC = max(mean(vFAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vFAll,2)==maxCC)-1)*tStep;
                      end
                      
                          
                      subplot(2,4,3+4*(vVar-1));
                      hold on;
                      scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),maxCC,...
                          40,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color, 'filled');
                      if colorID == 1
                        scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),maxCC,...
                              40,[0.5 0 0.5]);
                      else
                          scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),maxCC,...
                              40,[0 0.5 0]);
                      end
        %                   alpha(0.5);
                      ylim([-0.25 1]);
                      set(gca,'XTick',[1:4],'XTickLabels',{'dark','CL','CW','CCW'});
                      line([0 5],[0 0],'Color','k','LineStyle','--');

                      subplot(2,4,4+4*(vVar-1));
                      hold on;
                      scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),CClag,...
                          40,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color, 'filled');
                      if colorID == 1
                        scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),CClag,...
                              40,[0.5 0 0.5]);
                      else
                          scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),CClag,...
                              40,[0 0.5 0]);
                      end
                      ylim([-2.5 0]);
                      ylabel('peak time');
                      set(gca,'XTick',[1:4],'XTickLabels',{'dark','CL','CW','CCW'});
                  end
                  
              elseif contains(allAct{condID}.name,'PB')
                  if flyID == PBExFly
                    subplot(3,6,colorID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRPosRAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRPosLAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    ax = gca;
                    if colorID == 1
                        ylabel('vR > 0 autocorrelation');
                        text(-tStep*numPts/2-1,1.15,strcat(allAct{condID}.name,'-Summary'));
                        title('red channel');
                        ax.XColor = [0.5 0 0.5];
                        ax.YColor = [0.5 0 0.5];
                    else
                        title('green channel');
                        ax.XColor = [0 0.5 0];
                        ax.YColor = [0 0.5 0];
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 1], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 1]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                        if colorID == 1
                            legend({'dark (R)','dark (L)','CL (R)','CL(L)','CW (R)','CW (L)','CCW (R)','CCW (L)'});
                            legend('boxoff');
                        end
                    end
                    
                    subplot(3,6,6+colorID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRNegRAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRNegLAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    ax = gca;
                    if colorID == 1
                        ylabel('vR < 0 autocorrelation');
                        ax.XColor = [0.5 0 0.5];
                        ax.YColor = [0.5 0 0.5];
                    else
                        ax.XColor = [0 0.5 0.5];
                        ax.YColor = [0 0.5 0.5];
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 1], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 1]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                    end

                    subplot(3,6,12+colorID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFRAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFLAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                    line([0 0],[-0.25 1], 'Color', 'k','LineStyle','--');
                    ylim([-0.25 1]);
                    xlim([-tStep*numPts/2 tStep*numPts/2]);
                    ax = gca;
                    if colorID == 1
                        ylabel('vF autocorrelation');
                        ax.XColor = [0.5 0 0.5];
                        ax.YColor = [0.5 0 0.5];
                    else
                        ax.XColor = [0 0.5 0];
                        ax.YColor = [0 0.5 0];
                    end
                    xlabel('time delay (s)');
                  end
                  
                  for vVar = 1:6
                      if vVar == 1 
                          maxCC = max(mean(vRPosLAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vRPosLAll,2)==maxCC)-1)*tStep;
                      elseif vVar == 2
                          maxCC = max(mean(vRNegLAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vRNegLAll,2)==maxCC)-1)*tStep;
                      elseif vVar == 3
                          maxCC = max(mean(vFLAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vFLAll,2)==maxCC)-1)*tStep;
                      elseif vVar == 4
                          maxCC = max(mean(vRPosRAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vRPosRAll,2)==maxCC)-1)*tStep;
                      elseif vVar == 5
                          maxCC = max(mean(vRNegRAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vRNegRAll,2)==maxCC)-1)*tStep;
                      elseif vVar == 6
                          maxCC = max(mean(vFRAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vFRAll,2)==maxCC)-1)*tStep;
                      end
                      
                          
                      subplot(3,6,3+6*mod(vVar-1,3)+floor((vVar-1)/3));
                      hold on;
                      scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),maxCC,...
                          30,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color, 'filled');
                      if colorID == 1
                        scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),maxCC,...
                              30,[0.5 0 0.5]);
                      else
                          scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),maxCC,...
                              30,[0 0.5 0]);
                      end
        %                   alpha(0.5);
                      ylim([-0.25 1]);
                      set(gca,'XTick',[1:4],'XTickLabels',{'dark','CL','CW','CCW'});
                      line([0 5],[0 0],'Color','k','LineStyle','--');
                      if vVar < 4
                          title('left PB');
                      else
                          title('right PB');
                      end

                      subplot(3,6,5+6*mod(vVar-1,3)+floor((vVar-1)/3));
                      hold on;
                      scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),CClag,...
                          30,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color, 'filled');
                      if colorID == 1
                        scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),CClag,...
                              30,[0.5 0 0.5]);
                      else
                          scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),CClag,...
                              30,[0 0.5 0]);
                      end
                      ylim([-2.5 0]);
                      if vVar < 4
                          title('left PB');
                          ylabel('peak time');
                      else
                          title('right PB');
                      end
                      set(gca,'XTick',[1:4],'XTickLabels',{'dark','CL','CW','CCW'});
                  end
              end
          end
      end
   end
   reg.Renderer='Painters';
   set(reg,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
   print(reg,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\',...
   allAct{condID}.name,'\',allAct{condID}.name,'_CCSummary'),'-dpdf');
end

%% In the PB, look at the auto correlations for right and left turns for the PB only

% Look at the EB and PB
for condID = 2
   reg = figure('units','normalized','outerposition',[0 0 1 1]);
   
   numFlies = length(allAct{condID}.fly);
   % Step through the flies
   for flyID = 1:cond{condID}.numFlies
      % Sort the data by period
      for periodID = 2:4

          % Sort the data by color
          for colorID = 1
              if contains(allAct{condID}.name,'EB')
                  vRAll = [];
                  vFAll = [];
              elseif contains(allAct{condID}.name,'PB')
                  vRPosLAll = [];
                  vRNegLAll = [];
                  vFLAll = [];
                  vRPosRAll = [];
                  vRNegRAll = [];
                  vFRAll = [];
              end
              % Sort the data across trials
              for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

                  % Group across trials
                  if contains(allAct{condID}.name,'EB')
                      vRAll = horzcat(vRAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vRAutoC);
                      vFAll = horzcat(vFAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vFAutoC);
                  elseif contains(allAct{condID}.name,'PB') 
                      vRPosLAll = horzcat(vRPosLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRPosAutoC);
                      vRPosRAll = horzcat(vRPosLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRPosAutoC);

                      vRNegLAll = horzcat(vRNegLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRNegAutoC);
                      vRNegRAll = horzcat(vRNegLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRNegAutoC);

                      vFLAll = horzcat(vFLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vFAutoC);
                      vFRAll = horzcat(vFRAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vFAutoC);
                  end
              end
              if contains(allAct{condID}.name,'EB')

              elseif contains(allAct{condID}.name,'PB')
                    subplot(3,numFlies,flyID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRPosRAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRPosLAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    if colorID == 1
                        ylabel('vR > 0 autocorrelation');
                        if flyID == 1
                            text(-tStep*numPts/2,1.25,allAct{condID}.name);
                            line([0.4 0.6],[0.45 0.45],'Color','b');
                            text(0.7,0.45,'CL');
                            line([0.4 0.6],[0.4 0.4],'Color',[0.5 0 1]);
                            text(0.7,0.4,'CW');
                            line([0.4 0.6],[0.35 0.35],'Color',[0 0.5 1]);
                            text(0.7,0.35,'CCW');
                            line([0.4 0.6],[0.25 0.25],'Color','k');
                            text(0.7,0.25,'L PB');
                            line([0.4 0.6],[0.2 0.2],'Color','k','LineStyle','--');
                            text(0.7,0.2,'R PB');
                        end
                        title(strcat('fly #',num2str(flyID),'-PEG'));
                    else
                        title(strcat('fly #',num2str(flyID),'-green'));
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 0.5]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                    end
                    
                    subplot(3,numFlies,flyID+numFlies);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRNegRAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRNegLAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    if colorID == 1
                        ylabel('vR < 0 autocorrelation');
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 0.5]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                    end

                    subplot(3,numFlies,flyID+2*numFlies);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFRAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFLAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                    line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                    ylim([-0.25 0.5]);
                    xlim([-tStep*numPts/2 tStep*numPts/2]);
                    if colorID == 1 & flyID == 1
                        ylabel('vF autocorrelation');
                    end
                    xlabel('time delay (s)');
              end
          end
      end
   end                
end

set(reg,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(reg,'C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\PB\PEGPBAutoCC','-dpdf');

%% In the PB, look at the peak max vs. vRot

vRThresh = pi/8;

sgolayOrder = 3;
sgolayFrames = 11;

% Look at the PB condition specifically
for condID = 1:length(cond)
    if strcmp(allAct{condID}.name,'PB');
        for colorID = 1:2
            actR = figure('units','normalized','outerposition',[0 0 1 1]);
            for flyID = 1:length(allAct{condID}.fly)
                for periodID = 1:4
                    % Create arrays to hold max velocities and bump amplitudes
                    vRPAll = [];
                    vRNAll = [];
                    maxActLAllP = [];
                    maxActLAllN = [];
                    maxActRAllP = [];
                    maxActRAllN = [];
                    for trialID = 1:length(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial)
                        vR = sgolayfilt(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vR,...
                            sgolayOrder,sgolayFrames);
                        vRPos = vR(find(vR > vRThresh));
                        vRPAll = vertcat(vRPAll,vRPos);
                        vRNeg = vR(find(vR < -vRThresh));
                        vRNAll = vertcat(vRNAll,vRNeg);
                        
                        maxActL = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.act)-1;
                        maxActLAllP = vertcat(maxActLAllP,maxActL(find(vR > vRThresh))'); 
                        maxActLAllN = vertcat(maxActLAllN,maxActL(find(vR < -vRThresh))'); 
                        maxActR = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.act)-1;
                        maxActRAllP = vertcat(maxActRAllP,maxActR(find(vR > vRThresh))');
                        maxActRAllN = vertcat(maxActRAllN,maxActR(find(vR < -vRThresh))');
                        
                        subplot(2,numFlies,flyID);
                        hold on;
                        scatter(vR(find(abs(vR)>vRThresh)),maxActL(find(abs(vR)>vRThresh)),20,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,'filled');
                        alpha(0.1);
                        xlim([-pi pi]);
                        ylim([0 3-colorID]);
                        
                        subplot(2,numFlies,flyID+numFlies);
                        hold on;
                        scatter(vR(find(abs(vR)>vRThresh)),maxActR(find(abs(vR)>vRThresh)),20,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,'filled');
                        alpha(0.1);
                        xlim([-pi pi]);
                        ylim([0 3-colorID]);
                    end
                    [posFitR,posErrorR] = polyfit(vRPAll,maxActRAllP,1);
                    [negFitR,negErrorR] = polyfit(vRNAll,maxActRAllN,1);
                    [posFitL,posErrorL] = polyfit(vRPAll,maxActLAllP,1);
                    [negFitL,negErrorL] = polyfit(vRNAll,maxActLAllN,1);
                    
                    posV = linspace(vRThresh,pi,21);
                    negV = linspace(-pi,-vRThresh,21);
                    
                    yPosR = polyval(posFitR,posV);
                    yNegR = polyval(negFitR,negV);
                    yPosL = polyval(posFitL,posV);
                    yNegL = polyval(negFitL,negV);
                     
                    subplot(2,numFlies,flyID);
                    plot(posV,yPosL,'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(negV,yNegL,'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    
                    subplot(2,numFlies,flyID+numFlies);
                    plot(posV,yPosR,'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(negV,yNegR,'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                end
                if flyID == 1
                    subplot(2,numFlies,flyID)
                    text(-pi,3.1-colorID,'left PB');
                    ylabel('max DF/F');
                    subplot(2,numFlies,numFlies+flyID)
                    text(-pi,3.1-colorID,'right PB');
                    ylabel('max DF/F');
                end
                subplot(2,numFlies,flyID);
                title(strcat('fly #',num2str(flyID)));
                subplot(2,numFlies,numFlies+flyID)
                xlabel('vR (rad/s)');
            end
        end
    end
end

%% In the PB, look at the mean peak max vs. vRot - binned

% Specify rotational velocity bin edges and centers
vREdges = linspace(-pi,pi,9);
vRCents = vREdges(1:end-1)+0.5*mean(diff(vREdges));

% Specify Savtisky-Golay filter parameters for smoothing the velocity
sgolayOrder = 3;
sgolayFrames = 11;

% Look at the PB condition specifically
for condID = 1:length(cond)
    if strcmp(allAct{condID}.name,'PB')
        actMeanAll = zeros(2,4,2,length(vRCents));
        actSDAll = zeros(2,4,2,length(vRCents));
        for colorID = 1:2
            actBin = figure('units','normalized','outerposition',[0 0 1 1]);
            numFlies = length(allAct{condID}.fly);
            for flyID = 1:numFlies
                for periodID = 1:4
                    % Create arrays to hold max velocities and bump amplitudes
                    vRAll = [];
                    maxActLAll = [];
                    maxActRAll = [];
                    for trialID = 1:length(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial)
                        vR = sgolayfilt(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vR,...
                            sgolayOrder,sgolayFrames);
                        vRAll = vertcat(vRAll,vR);
                        
                        maxActL = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.act)-1;
                        maxActLAll = vertcat(maxActLAll,maxActL'); 
                        maxActR = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.act)-1;
                        maxActRAll = vertcat(maxActRAll,maxActR');
                    end
                    binIDs = discretize(vRAll,vREdges);
                    for binNow = 1:length(vREdges)-1
                        meanAct = mean(maxActLAll(find(binIDs == binNow)));
                        actMeanAll(colorID,periodID,1,binNow) = meanAct;
                        SDAct = std(maxActLAll(find(binIDs == binNow)));
                        actSDAll(colorID,periodID,1,binNow) = SDAct;
                    end
                    
                    for binNow = 1:length(vREdges)-1
                        meanAct = mean(maxActRAll(find(binIDs == binNow)));
                        actMeanAll(colorID,periodID,2,binNow) = meanAct;
                        SDAct = std(maxActRAll(find(binIDs == binNow)));
                        actSDAll(colorID,periodID,2,binNow) = SDAct;
                    end
                end
                
                subplot(3,numFlies,flyID);
                hold on;
                title(strcat('fly #',num2str(flyID)));
                for periodID = 1:4
                    meanNow = squeeze(actMeanAll(colorID,periodID,1,:));
                    SDNow = squeeze(actSDAll(colorID,periodID,1,:));
                    plot(vRCents,meanNow,...
                        'color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    p = patch('XData',[vRCents fliplr(vRCents)],...
                        'YData',[meanNow'+SDNow' fliplr(meanNow'-SDNow')],...
                        'FaceColor',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                        'FaceAlpha',0.1,...
                        'EdgeColor','none');
                    uistack(p,'bottom');
                    
                end
                xlim([vREdges(1) vREdges(end)]);
                ylim([0.25 1.25-0.5*(colorID-1)]);
                line([0 0],[0.25 1.25-0.5*(colorID-1)],'Color','k','LineStyle','--');
                
                subplot(3,numFlies,numFlies+flyID);
                hold on;
                for periodID = 1:4
                    meanNow = squeeze(actMeanAll(colorID,periodID,2,:));
                    SDNow = squeeze(actSDAll(colorID,periodID,2,:));
                    plot(vRCents,meanNow,...
                        'color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    p = patch('XData',[vRCents fliplr(vRCents)],...
                        'YData',[meanNow'+SDNow' fliplr(meanNow'-SDNow')],...
                        'FaceColor',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                        'FaceAlpha',0.1,...
                        'EdgeColor','none');
                    uistack(p,'bottom');
                end
                xlim([vREdges(1) vREdges(end)]);
                ylim([0.25 1.25-0.5*(colorID-1)]);
                line([0 0],[0.25 1.25-0.5*(colorID-1)],'Color','k','LineStyle','--');
                
                subplot(3,numFlies,2*numFlies+flyID);
                hold on;
                for periodID = 1:4
                    meanNowL = squeeze(actMeanAll(colorID,periodID,1,:));
                    meanNowR = squeeze(actMeanAll(colorID,periodID,2,:));
                    plot(vRCents,meanNowR-meanNowL,...
                        'color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                end
                xlabel('vR (rad/s)');
                xlim([vREdges(1) vREdges(end)]);
                ylim([-0.15 0.15]);
                line([0 0],[-0.15 0.15],'Color','k','LineStyle','--');
                line([vREdges(1) vREdges(end)],[0 0],'Color','k','LineStyle','--');
                
                if flyID == 1
                    subplot(3,numFlies,flyID)
                    text(vREdges(1),1.3-0.5*(colorID-1),'left PB');
                    if colorID == 1
                        text(vREdges(1)-pi/4,1.4-0.5*(colorID-1),'red channel');
                    else
                        text(vREdges(1)-pi/4,1.4-0.5*(colorID-1),'green channel');
                    end
                    ylabel('max DF/F');
                    subplot(3,numFlies,numFlies+flyID)
                    text(vREdges(1),1.3-0.5*(colorID-1),'right PB');
                    ylabel('max DF/F');
                    subplot(3,numFlies,2*numFlies+flyID)
                    text(vREdges(1),0.175,'left - right');
                    ylabel('DF/F');
                end
            end
            set(actBin,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
            if colorID == 1
                print(actBin,'C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\PB\PBActVSVR_Red','-dpdf');
            else
                print(actBin,'C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\PB\PBActVSVR_Green','-dpdf');
            end
        end
    end
end

%% In the PB, look at the mean peak max vs. vF - binned

% Specify rotational velocity bin edges and centers
vFEdges = linspace(0,1,9);
vFCents = vFEdges(1:end-1)+0.5*mean(diff(vFEdges));

% Specify Savtisky-Golay filter parameters for smoothing the velocity
sgolayOrder = 3;
sgolayFrames = 11;

% Look at the PB condition specifically
for condID = 1:length(cond)
    if strcmp(allAct{condID}.name,'PB')
        actMeanAll = zeros(2,4,2,length(vFCents));
        actSDAll = zeros(2,4,2,length(vFCents));
        for colorID = 1:2
            actBin = figure('units','normalized','outerposition',[0 0 1 1]);
            numFlies = length(allAct{condID}.fly);
            for flyID = 1:numFlies
                for periodID = 1:4
                    % Create arrays to hold max velocities and bump amplitudes
                    vFAll = [];
                    maxActLAll = [];
                    maxActRAll = [];
                    for trialID = 1:length(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial)
                        vF = sgolayfilt(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vF,...
                            sgolayOrder,sgolayFrames);
                        vFAll = vertcat(vFAll,vF);
                        
                        maxActL = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.act)-1;
                        maxActLAll = vertcat(maxActLAll,maxActL'); 
                        maxActR = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.act)-1;
                        maxActRAll = vertcat(maxActRAll,maxActR');
                    end
                    binIDs = discretize(vFAll,vFEdges);
                    for binNow = 1:length(vFEdges)-1
                        meanAct = mean(maxActLAll(find(binIDs == binNow)));
                        actMeanAll(colorID,periodID,1,binNow) = meanAct;
                        SDAct = std(maxActLAll(find(binIDs == binNow)));
                        actSDAll(colorID,periodID,1,binNow) = SDAct;
                    end
                    
                    for binNow = 1:length(vFEdges)-1
                        meanAct = mean(maxActRAll(find(binIDs == binNow)));
                        actMeanAll(colorID,periodID,2,binNow) = meanAct;
                        SDAct = std(maxActRAll(find(binIDs == binNow)));
                        actSDAll(colorID,periodID,2,binNow) = SDAct;
                    end
                end
                
                subplot(3,numFlies,flyID);
                hold on;
                title(strcat('fly #',num2str(flyID)));
                for periodID = 1:4
                    meanNow = squeeze(actMeanAll(colorID,periodID,1,:));
                    SDNow = squeeze(actSDAll(colorID,periodID,1,:));
                    plot(vFCents,meanNow,...
                        'color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    p = patch('XData',[vFCents fliplr(vFCents)],...
                        'YData',[meanNow'+SDNow' fliplr(meanNow'-SDNow')],...
                        'FaceColor',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                        'FaceAlpha',0.1,...
                        'EdgeColor','none');
                    uistack(p,'bottom');
                    
                end
                xlim([vFEdges(1) vFEdges(end)]);
                ylim([0.25 1.25-0.5*(colorID-1)]);
                line([0 0],[0.25 1.25-0.5*(colorID-1)],'Color','k','LineStyle','--');
                
                subplot(3,numFlies,numFlies+flyID);
                hold on;
                for periodID = 1:4
                    meanNow = squeeze(actMeanAll(colorID,periodID,2,:));
                    SDNow = squeeze(actSDAll(colorID,periodID,2,:));
                    plot(vFCents,meanNow,...
                        'color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    p = patch('XData',[vFCents fliplr(vFCents)],...
                        'YData',[meanNow'+SDNow' fliplr(meanNow'-SDNow')],...
                        'FaceColor',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                        'FaceAlpha',0.1,...
                        'EdgeColor','none');
                    uistack(p,'bottom');
                end
                xlim([vFEdges(1) vFEdges(end)]);
                ylim([0.25 1.25-0.5*(colorID-1)]);
                line([0 0],[0.25 1.25-0.5*(colorID-1)],'Color','k','LineStyle','--');
                
                subplot(3,numFlies,2*numFlies+flyID);
                hold on;
                for periodID = 1:4
                    meanNowL = squeeze(actMeanAll(colorID,periodID,1,:));
                    meanNowR = squeeze(actMeanAll(colorID,periodID,2,:));
                    plot(vFCents,meanNowR-meanNowL,...
                        'color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                end
                xlabel('vF (cm/s)');
                xlim([vFEdges(1) vFEdges(end)]);
                ylim([-0.15 0.15]);
                line([0 0],[-0.15 0.15],'Color','k','LineStyle','--');
                line([vFEdges(1) vFEdges(end)],[0 0],'Color','k','LineStyle','--');
                
                if flyID == 1
                    subplot(3,numFlies,flyID)
                    text(vFEdges(1),1.3-0.5*(colorID-1),'left PB');
                    if colorID == 1
                        text(vFEdges(1)-pi/4,1.4-0.5*(colorID-1),'red channel');
                    else
                        text(vFEdges(1)-pi/4,1.4-0.5*(colorID-1),'green channel');
                    end
                    ylabel('max DF/F');
                    subplot(3,numFlies,numFlies+flyID)
                    text(vFEdges(1),1.3-0.5*(colorID-1),'right PB');
                    ylabel('max DF/F');
                    subplot(3,numFlies,2*numFlies+flyID)
                    text(vFEdges(1),0.175,'left - right');
                    ylabel('DF/F');
                end
            end
            set(actBin,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
            if colorID == 1
                print(actBin,'C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\PB\PBActVSVF_Red','-dpdf');
            else
                print(actBin,'C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\PB\PBActVSVF_Green','-dpdf');
            end
        end
    end
end

%% In the EB, subtract the PEG bump from the PEN2s and plot

% Look at the EB data
reg = figure('units','normalized','outerposition',[0 0 1 1]);

% Set a threshold for when to find the PEG max
actThresh = 1.25;

% Find the angles around the EB
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs+0.5*mean(diff(angs));

% Step through each fly and trial
for condID = 1:length(cond)
   if contains(cond{condID}.name,'EB')
       % Step through the flies
       for flyID = 1:cond{condID}.numFlies

          % Sort the data across trials
          for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

             % Extract the behavioral parameters
             tPts = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.OffsetRotMatch(:,1);
             heading = pi/180*cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.PosRotMatch;

             % Find the max PEG ROI values and set the PEN values around there to
             % 0
             maxPEG = max(cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax);
             PEN2SubPlt = cond{condID}.allFlyData{flyID}.All{trialID}.GROIaveMax-1;
             for maxFnd = 1:length(maxPEG)
                 if maxPEG(maxFnd) > actThresh
                     pkPt = find(maxPEG(maxFnd) == cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax(:,maxFnd));
                     PEN2SubPlt(pkPt,maxFnd) = -0.1;
                     if pkPt > 1
                         PEN2SubPlt(pkPt-1,maxFnd) = -0.1;
                     else
                         PEN2SubPlt(16,maxFnd) = -0.1;
                     end
                     if pkPt < 16
                         PEN2SubPlt(pkPt+1,maxFnd) = -0.1;
                     else
                         PEN2SubPlt(1,maxFnd) = -0.1;
                     end
                 end
             end
             
             % Plot the profile
             subplot(5,cond{condID}.numFlies,flyID+cond{condID}.numFlies*(trialID-1));
             hold on;
             imagesc(tPts,angs,PEN2SubPlt);
%              plot(tPts,-heading,'k','LineStyle','--');
             colormap(brewermap(64, 'Blues'));
             caxis([-0.1 1.5]);
             colorbar;
             ylim([-pi pi]);
             xlim([tPts(1) tPts(end)]);
             if flyID == 1
                 ylabel('EB position (rad)');
             end
             if trialID == 1
                 title(strcat('fly #',num2str(flyID)));
             end
             if trialID == 5
                 xlabel('time (s)');
             end
          end
       end
       set(reg,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
        print(reg,'C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\EB\EB_PEN2-PEGBump','-dpdf');
   end
end

%% In the EB, register the peaks with respect to one another and compare

% Look at the EB data
reg = figure('units','normalized','outerposition',[0 0 1 1]);

% Set a threshold for when to find the PEG max
actThreshR = 1.25;
actThreshG = 1.5;

% Set a velocity threshold
vRThresh = 0;

% Find the angles around the EB
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs+0.5*mean(diff(angs));

% Step through each fly and trial
for condID = 1:length(cond)
   if contains(cond{condID}.name,'EB')
       % Step through the flies
       for flyID = 1:cond{condID}.numFlies

          % Sort the data across trials
          for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

             % Extract the behavioral parameters
             tPts = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.OffsetRotMatch(:,1);
             vR = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.vRot;
             vR = sgolayfilt(vR,3,11);

             % Find the max PEG ROI values and set the PEN values around there to
             % 0
             maxPEG = max(cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax);
             maxPEN2 = max(cond{condID}.allFlyData{flyID}.All{trialID}.GROIaveMax);
             
             PEGReg = zeros(size(cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax));
             PEN2Reg = zeros(size(cond{condID}.allFlyData{flyID}.All{trialID}.GROIaveMax));
             
             for maxFnd = 1:length(maxPEG)
                 if maxPEG(maxFnd) > actThreshR
                     XCR = cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax(:,maxFnd);
                     XCG = cond{condID}.allFlyData{flyID}.All{trialID}.GROIaveMax(:,maxFnd);
                     pkPt = find(maxPEG(maxFnd) == XCR);
                     PEN2Reg(:,maxFnd) = circshift(XCG-1,8-pkPt);
                 end
                 if maxPEN2(maxFnd) > actThreshG
                     XCR = cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax(:,maxFnd);
                     XCG = cond{condID}.allFlyData{flyID}.All{trialID}.GROIaveMax(:,maxFnd);
                     pkPt = find(maxPEN2(maxFnd) == XCG);
                     PEGReg(:,maxFnd) = circshift(XCR-1,8-pkPt);
                 end
             end
             
             vRHigh = find(abs(vR) > vRThresh);
             vRPos = find(vR > 0);
             vRNeg = find(vR < 0);
             
             % Plot the profile
             subplot(5,cond{condID}.numFlies,flyID+cond{condID}.numFlies*(trialID-1));
             hold on;
             plot(angs,(mean(PEGReg(:,vRHigh),2)-min(mean(PEGReg(:,vRHigh),2)))./...
                 (max(mean(PEGReg(:,vRHigh),2))-min(mean(PEGReg(:,vRHigh),2))),'r');
%              plot(angs,mean(PEGReg(:,vRPos),2)./max(mean(PEGReg(:,vRPos),2)),'r');
%              plot(angs,mean(PEGReg(:,vRNeg),2)./max(mean(PEGReg(:,vRNeg),2)),'r',...
%                  'LineWidth',2','LineStyle','--');
             plot(angs,(mean(PEN2Reg(:,vRHigh),2)-min(mean(PEN2Reg(:,vRHigh),2)))./...
                 (max(mean(PEN2Reg(:,vRHigh),2))-min(mean(PEN2Reg(:,vRHigh),2))),'g');
%              plot(angs,mean(PEN2Reg(:,vRPos),2)./max(mean(PEN2Reg(:,vRPos),2)),'g');
%              plot(angs,mean(PEN2Reg(:,vRNeg),2)./max(mean(PEN2Reg(:,vRNeg),2)),'g',...
%                  'LineWidth',2','LineStyle','--');
             xlim([-pi pi]);
          end
       end
   end
end

%% Look at the bump persistance

% Set the Savitsky Golay filter parameters
sgolayOrder = 3;
sgolayFrames = 11;

% Set the velocity thresholds that define standing bouts
vRThresh = 0.01*pi;
vFThresh = 0.01;

% Set a threshold for the length of a standing bout
minBoutL = 5;

% Step through the conditions
for condID = 1:length(cond)
    % Step through the flies
    for flyID = 1:cond{condID}.numFlies
        
        % Step through the colors
        for colorID = 1:2

            numPeriods = length(allAct{condID}.fly{flyID}.color{colorID}.period);
            
            % Sort the data by period
            for periodID = 1:numPeriods

                % Inititalize variables
                maxBoutL = 0; % Max bout length
                numBouts = 0; % Number of bouts

                % Sort the data across trials
                for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

                    % Find the time step
                    tPts = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.tPts;
                    tStep = mean(diff(tPts));

                    % Find when the animal is standing
                    vR = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vR;
                    vR = sgolayfilt(vR,sgolayOrder,sgolayFrames);
                    vF = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vF;
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

                    if ~isempty(boutStart)
                        maxBoutL = max(maxBoutL,max(boutStop-boutStart)+1);
                    end
                    numBouts = numBouts+length(boutStart);

                    % Add the standing bout info to the object
                    allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.boutStart = standingBouts(boutStart);
                    allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.boutStop = standingBouts(boutStop);

                end

                allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.maxBoutL = maxBoutL;
                allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.numBouts = numBouts;
            end
        end
    end
end

% Group the standing bouts across trials

% Make the figure
bumpPerEB = figure('units','normalized','outerposition',[0 0 1 1]);
bumpPerPB = figure('units','normalized','outerposition',[0 0 1 1]);

for condID = 1:length(cond)
    % Step through the flies
    for flyID = 1:cond{1}.numFlies

        % Step through the colors
        for colorID = 1:2

            numPeriods = length(allAct{condID}.fly{flyID}.color{colorID}.period);
            % Sort the data by period
            for periodID = 1:numPeriods

                % Pull out the relevant values
                maxBoutL = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.maxBoutL;
                numBouts = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.numBouts;

                standBout = 1;
                if contains(allAct{condID}.name,'EB')
                    standAct = zeros(numBouts,16,maxBoutL);
                    standActMean = zeros(16,maxBoutL);
                elseif contains(allAct{condID}.name,'PB')
                    standActL = zeros(numBouts,8,maxBoutL);
                    standActR = zeros(numBouts,8,maxBoutL);
                    standActMeanL = zeros(8,maxBoutL);
                    standActMeanR = zeros(8,maxBoutL);
                end

                % Sort the data across trials
                for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

                    if contains(allAct{condID}.name,'EB')
                        % Place all of the standing bout activity into the array
                        actNow = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.act-1;
                        boutStart = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.boutStart;
                        boutStop = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.boutStop;

                        if ~isempty(boutStart)
                            for bout = 1:length(boutStart)
                                actBout = actNow(:,boutStart(bout):boutStop(bout));
                                maxROI = find(actBout(:,1) == max(actBout(:,1)));
                                actBoutShift = circshift(actBout,8-maxROI,1);
                                standAct(standBout,:,1:(boutStop(bout)-boutStart(bout)+1)) = actBoutShift;
                                standBout = standBout + 1;
                            end
                        end
                    elseif contains(allAct{condID}.name,'PB')
                        % Place all of the standing bout activity into the array
                        actNowL = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.act-1;
                        actNowR = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.act-1;
                        boutStart = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.boutStart;
                        boutStop = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.boutStop;

                        if ~isempty(boutStart)
                            for bout = 1:length(boutStart)
                                actBoutL = actNowL(:,boutStart(bout):boutStop(bout));
                                actBoutR = actNowR(:,boutStart(bout):boutStop(bout));
                                maxROIL = find(actBoutL(:,1) == max(actBoutL(:,1)));
                                maxROIR = find(actBoutR(:,1) == max(actBoutR(:,1)));
                                actBoutShiftL = circshift(actBoutL,4-maxROIL,1);
                                actBoutShiftR = circshift(actBoutR,4-maxROIR,1);
                                standActL(standBout,:,1:(boutStop(bout)-boutStart(bout)+1)) = actBoutShiftL;
                                standActR(standBout,:,1:(boutStop(bout)-boutStart(bout)+1)) = actBoutShiftR;
                                standBout = standBout + 1;
                            end
                        end
                    end
                end

                % Find the mean activity profile
                for l = 1:maxBoutL
                    if contains(allAct{condID}.name,'EB')
                        allActNow = squeeze(standAct(:,:,l));
                        zeroPers = [];
                        for bout = 1:size(allActNow,1)
                            if max(allActNow(bout,:)) == 0
                                zeroPers = vertcat(zeroPers,bout);
                            end
                        end
                        allActNow(zeroPers,:) = [];
                        standActMean(:,l) = mean(allActNow,1);
                    elseif contains(allAct{condID}.name,'PB')
                        allActNowL = squeeze(standActL(:,:,l));
                        allActNowR = squeeze(standActR(:,:,l));
                        zeroPers = [];
                        for bout = 1:size(allActNowL,1)
                            if max(allActNowL(bout,:)) == 0
                                zeroPers = vertcat(zeroPers,bout);
                            end
                        end
                        allActNowL(zeroPers,:) = [];
                        allActNowR(zeroPers,:) = [];
                        standActMeanL(:,l) = mean(allActNowL,1);
                        standActMeanR(:,l) = mean(allActNowR,1);
                    end
                end

                % Plot the mean values
                angVals = linspace(-pi,pi,17);
                angVals(end) = [];
                angVals = angVals + 0.5*mean(diff(angVals));
                if contains(allAct{condID}.name,'EB')
                    allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.standAct = standAct;
                    allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.standActMean = standActMean;
                    
                    figure(bumpPerEB);
                    subplot(numPeriods,2*cond{condID}.numFlies,2*flyID+colorID-2+2*cond{condID}.numFlies*(periodID-1))
                    tVals = tStep*[0:size(standActMean,2)-1];
                    imagesc(tVals,angVals,standActMean);
                    if flyID == 1
                        ylabel('EB position (rad');
                    end
                    if colorID == 1
                        caxis([0 1]);
                    else
                        caxis([0 0.75]);
                    end
                elseif contains(allAct{condID}.name,'PB')
                    allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.standActL = standActL;
                    allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.standActR = standActR;
                    allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.standActMeanL = standActMeanL;
                    allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.standActMeanR = standActMeanR;
                    
                    figure(bumpPerPB);
                    subplot(numPeriods,4*cond{condID}.numFlies,4*flyID+colorID-4+4*cond{condID}.numFlies*(periodID-1))
                    tVals = tStep*[0:size(standActMeanL,2)-1];
                    imagesc(tVals,[1:8],standActMeanL);
                    if flyID == 1
                        ylabel('PB position (rad');
                    end
                    xlim([0 5]);
                    if colorID == 1
                        caxis([0 1]);
                    else
                        caxis([0 0.75]);
                    end
                    colormap(brewermap(64, 'Blues'));
                    if periodID == 4
                        xlabel('time (s)');
                    end
                    if periodID == 1
                       text(0,-1.25*pi,strcat('fly #',num2str(flyID))); 
                    end
                    title(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type);
                    
                    subplot(numPeriods,4*cond{condID}.numFlies,4*flyID+colorID-2+4*cond{condID}.numFlies*(periodID-1))
                    tVals = tStep*[0:size(standActMeanR,2)-1];
                    imagesc(tVals,[1:8],standActMeanR);
                    if colorID == 1
                        caxis([0 1]);
                    else
                        caxis([0 0.75]);
                    end
                end
                xlim([0 5]);
                colormap(brewermap(64, 'Blues'));
                if periodID == 4
                    xlabel('time (s)');
                end
                if periodID == 1
                   text(0,-1.25*pi,strcat('fly #',num2str(flyID))); 
                end
                title(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type);
            end
        end
    end
end

set(bumpPerEB,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpPerEB,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\EB\EB_Persistence'),'-dpdf');

set(bumpPerPB,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpPerPB,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\PB\PB_Persistence'),'-dpdf');

%% Look at the bump persistance - create nice summary figure

% Example flies and standing bouts
EBFlyEx = 2;
EBBout = 3;

PBFlyEx = 5;
PBBout = 1;

% Make the figure
bumpPer = figure('units','normalized','outerposition',[0 0 1 1]);

% Step through the conditions
for condID = 1:length(cond)
    % Step through the flies
    for flyID = 1:cond{1}.numFlies

        % Step through the colors
        for colorID = 1:2

            % Sort the data by period
            for periodID = 1

                if contains(allAct{condID}.name,'EB')
                    
                    % Specify the EB angles
                    angVals = linspace(-pi,pi,17);
                    angVals(end) = [];
                    angVals = angVals + 0.5*mean(diff(angVals));
                    
                    % Pull out the activity
                    standAct = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.standAct;
                    allStandAct = squeeze(standAct(:,8,:));
                    
                    % Find the time step
                    tPts = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.tPts;
                    tStep = mean(diff(tPts));
                    tVals = tStep*[0:size(standAct,3)-1];
                    
                    % Normalize each trial for plotting and get rid of zero
                    % values from the array initialization
                    clear standActMean;
                    for l = 1:allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.maxBoutL
                        allActNow = squeeze(standAct(:,:,l));
                        zeroPers = [];
                        for bout = 1:size(allActNow,1)
                            if max(allActNow(bout,:)) == 0
                                zeroPers = vertcat(zeroPers,bout);
                            else
                                allActNow(bout,:) = allActNow(bout,:)./standAct(bout,8,1);
                            end
                        end
                        allActNow(zeroPers,:) = [];
                        standActMean(:,l) = mean(allActNow,1);
                    end
                    
                    % Plot the EB example
                    if flyID == EBFlyEx
                        subplot(4,4,4*colorID-3);
                        hold on;
                        imagesc(tVals,angVals,squeeze(standAct(EBBout,:,:)));
                        line([0 3],[angVals(8) angVals(8)],'Color','k','LineStyle','--');
                        caxis([0 1]);
                        xlim([0 3]);
                        ylim([-pi pi]);
                        ylabel('EB position (rad)');
                        if colorID == 1
                            title('red channel');
                        else
                            title('green channel');
                        end
                        colorbar;
                        
                        subplot(4,4,4*colorID-2);
                        plot(tVals,squeeze(standAct(EBBout,8,:)));
                        xlim([0 3]);
                        ylim([0 1]);
                        ylabel('DF/F');
                        
                        subplot(4,4,4*colorID-1);
                        hold on;
                        for bout = 1:size(standAct,1)
                            allStandAct(bout,:) = allStandAct(bout,:)./allStandAct(bout,1);
                        end
                        plot(tVals,allStandAct);
                        plot(tVals,squeeze(standActMean(8,:)),'k','LineWidth',2);
                        f = fit(tVals',squeeze(standActMean(8,:))','exp1');
                        plot(tVals,f.a*exp(f.b*tVals),'Color','r');
                        text(2,1.25,strcat('tau =',num2str(-1/f.b),' s'));
                        xlim([0 3]);
                        ylabel('normalized F');
                        ylim([0 2]);
                        
                        colormap(brewermap(64, 'Blues'));
                    end
                    
                    % Calculate and plot the decay time constant
                    subplot(4,4,[4 8]);
                    hold on;
                    f = fit(tVals',squeeze(standActMean(8,:))','exp1');
                    scatter(colorID+0.05*(flyID-3),-1/f.b,40,[2-colorID colorID-1 2-colorID],'filled');
                    xlim([0 3]);
                    set(gca,'XTick',[1 2],'XTickLabels',{'red','green'});
                    ylim([0 20]);
                    ylabel('tau (s)');
                  
                    
                elseif contains(allAct{condID}.name,'PB')
                    
                    % Pull out the activity
                    standActL = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.standActL;
                    standActR = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.standActR;
                    allStandActL = squeeze(standActL(:,4,:));
                    allStandActR = squeeze(standActR(:,4,:));
                    
                    % Find the time step
                    tPts = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.tPts;
                    tStep = mean(diff(tPts));
                    tVals = tStep*[0:size(standActL,3)-1];
                    
                    % Normalize the data for plotting and remove 0 values
                    clear standActMeanL;
                    clear standActMeanR;
                    for l = 1:allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.maxBoutL
                        allActNowL = squeeze(standActL(:,:,l));
                        allActNowR = squeeze(standActR(:,:,l));
                        zeroPers = [];
                        for bout = 1:size(allActNowL,1)
                            if max(allActNowL(bout,:)) == 0
                                zeroPers = vertcat(zeroPers,bout);
                            else
                                allActNowL(bout,:) = allActNowL(bout,:)./standActL(bout,4,1);
                                allActNowR(bout,:) = allActNowR(bout,:)./standActR(bout,4,1);
                            end
                        end
                        allActNowL(zeroPers,:) = [];
                        allActNowR(zeroPers,:) = [];
                        standActMeanL(:,l) = mean(allActNowL,1);
                        standActMeanR(:,l) = mean(allActNowR,1);
                    end
                    
                    % Plot the PB example
                    if flyID == PBFlyEx
                        subplot(4,4,8+4*colorID-3);
                        hold on;
                        imagesc(tVals,[1:8],squeeze(standActL(PBBout,:,:)));
                        line([0 3],[4 4],'Color','k','LineStyle','--');
                        caxis([0 1]);
                        xlim([0 3]);
                        ylim([0 9]);
                        ylabel('PB position (glom)');
                        if colorID == 1
                            title('red channel');
                        else
                            title('green channel');
                            xlabel('standing time (s)');
                        end
                        colorbar;
                        
                        subplot(4,4,8+4*colorID-2);
                        plot(tVals,squeeze(standActL(PBBout,4,:)));
                        xlim([0 3]);
                        ylim([0 1]);
                        ylabel('DF/F');
                        if colorID == 2
                            xlabel('standing time (s)');
                        end
                        
                        subplot(4,4,8+4*colorID-1);
                        hold on;
                        for bout = 1:size(standActL,1)
                            allStandActL(bout,:) = allStandActL(bout,:)./allStandActL(bout,1);
                        end
                        plot(tVals,allStandActL);
                        plot(tVals,squeeze(standActMeanL(4,:)),'k','LineWidth',2);
                        f = fit(tVals',squeeze(standActMeanL(4,:))','exp1');
                        plot(tVals,f.a*exp(f.b*tVals),'Color','r');
                        text(2,1.25,strcat('tau =',num2str(-1/f.b),' s'));
                        xlim([0 3]);
                        ylabel('normalized F');
                        if colorID == 2
                            xlabel('standing time (s)');
                        end
                        ylim([0 2]);
                        
                        colormap(brewermap(64, 'Blues'));
                    end
                    
                    % Find and plot the exponential coefficient
                    subplot(4,4,[12 16]);
                    hold on;
                    
                    f1 = fit(tVals',squeeze(standActMeanL(4,:))','exp1');
                    scatter(colorID-0.125,-f1.b,40,[2-colorID colorID-1 2-colorID],'filled');
                    f2 = fit(tVals',squeeze(standActMeanR(4,:))','exp1');
                    scatter(colorID+0.125,-f2.b,40,[2-colorID colorID-1 2-colorID],'filled');
                    line([colorID-0.125 colorID+0.125],[-f1.b -f2.b],'color',[0.9 0.9 0.9]);
                    xlim([0 3]);
                    set(gca,'XTick',[1 2],'XTickLabels',{'red','green'});
%                     ylim([-0.4 0.4]);
                    ylabel('1/tau Hz');
                    line([0 3],[0 0],'color','k','LineStyle','--');
                    
                end
            end
        end
    end
end

set(bumpPer,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(bumpPer,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\PersistenceSummary'),'-dpdf');
