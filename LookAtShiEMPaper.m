%% Clear out old data
clear;
clc;

%60 s dark, 60 s closed, jump, 30s, jump, 30s, 30s CW, 30s CCW
%all data contained together but can sort by period

dir = '/Users/loaner/Documents/imaging/Data_Dan/shi/empty/'
shi = 'empty'

%dir = '/Users/loaner/Documents/imaging/Data_Dan/shi/EPG/'
%shi = 'EPG'

%dir = '/Users/loaner/Documents/imaging/Data_Dan/shi/GE/'
%shi = 'GE'

%dir = '/Users/loaner/Documents/imaging/Data_Dan/shi/PEG/'
%shi = 'PEG'

%dir = '/Users/loaner/Documents/imaging/Data_Dan/shi/PEN2/'
%shi = 'PEN2'

%dir = '/Users/loaner/Documents/imaging/Data_Dan/shi/PEN1/'
%shi = 'PEN1'


try
    from_file = load(strcat(dir, 'cond'), 'cond');
    cond=from_file.cond;    
catch
    cond = FlyDatLoad(1);
    save(strcat(dir, 'cond'), 'cond');
end



%cond = FlyDatLoad(1);

%% Look at the processed data across the various visual conditons for each fly

% Step across conditions
for condID = 1:length(cond)
    
    % Step across flies
    for flyID = 1:cond{condID}.numFlies
        
        % Make a figure
        actPlotCL = figure('units','normalized','outerposition',[0 0 1 1]);
        actPlotOL = figure('units','normalized','outerposition',[0 0 1 1]);
        actPlotDirCL = figure('units','normalized','outerposition',[0 0 1 1]);
        actPlotvRot = figure('units','normalized','outerposition',[0 0 1 1]);    
        
        % Step across hot and cold conditions
        for hc = 1:2
            if hc == 1
                trialType = 'All';
                TCol = 'b';
            else
                trialType = 'All_30C';
                TCol = 'r';
            end
            
            for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialType))
            
                % Sort out the closed loop and dark periods
                [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch);
                darkEnd = find(diff(darkPer)>1);
                darkPer(darkEnd:end) = [];
                
                % Pull out the relevant data
                tPts = cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.OffsetRotMatch(:,1);
                tPts = tPts - tPts(1);
                heading = pi/180*cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.PosRotMatch;
                stripePos = cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.OffsetRotMatch(:,2);
                DF = cond{condID}.allFlyData{flyID}.(trialType){trialID}.ROIaveMax-1;
                vRot = cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.vRot;
                
                %% Get directions and rate of change
               
                
                sgolayOrder = 3;
                sgolayFrames = 7;
                s = size(DF);
                newDF = [];
                for ind = 1:s(1);
                    g = DF(ind,:);
                    newDF(ind,:) = sgolayfilt(g,sgolayOrder,sgolayFrames);
                end
                vRot = sgolayfilt(vRot,sgolayOrder,sgolayFrames);
                    
                mags = [];
                dirs = [];
                for ind = 1:length(tPts);
                    vec = DF(:, ind);
                    [di, mi] = getVecSum( vec );
                    mags = [mags, mi/sum(abs(vec))];
                    dirs = [dirs di];
                end
                    
                    
                dirs = dirs(mags > 0.25); %only consider timepoints where directionality is reasonably strong
                tPtsVel = tPts(mags > 0.25);
                dirsPlt = dirs-pi; %consider range from -pi to pi
                vDF = transpose( getRates(dirsPlt, tPtsVel, 2*pi) ); %get rate of change of direction of fluorescent intensity
                
                %% Get heading Data
                % Convert jumps in heading or stripe position to NaNs for
                % plotting
                headingPlt = heading;
                stripePosPlt = stripePos;
                for step = 2:length(heading)
                    if abs(heading(step)-heading(step-1)) > pi
                        headingPlt(step-1) = NaN;
                    end
                    if abs(stripePos(step)-stripePos(step-1)) > pi
                        stripePosPlt(step-1) = NaN;
                    end
                end
                for step = 2:length(dirs)
                    if abs(dirs(step)-dirs(step-1)) > pi
                        dirsPlt(step-1) = NaN;
                    end
                end
                
                plotPer = vertcat(darkPer,CLPer);
                stripeJumps = find(diff(cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.Trans) ~= 0);
                stripeJumps(end) = [];
                
                tPtsDir = plotPer(mags(plotPer) > 0.25); %for directionality
                plotPerDir = [1 : length(tPtsDir)]; %for directionality
                
                
                %% Plot CL
                figure(actPlotCL);
                % Plot the dark periods
                subplot(5,1,trialID+2*(hc-1));
                hold on;
                imagesc(tPts(plotPer),linspace(-pi,pi,size(DF,1)),DF(:,plotPer));
                plot(tPts(darkPer),headingPlt(darkPer),'k');
                plot(tPts(CLPer),stripePosPlt(CLPer),'b');
                for jump = 1:length(stripeJumps);
                    line([tPts(stripeJumps(jump)) tPts(stripeJumps(jump))],...
                        [-pi pi],'Color','m');
                end
                axis tight;
                caxis([0 1]);
                if hc == 1 & trialID == 1
                    text(0,1.5*pi,strcat(cond{condID}.name,' - fly #',num2str(flyID)),'FontSize',14);
                    text(tPts(stripeJumps(1)),1.25*pi,'jump #1','FontSize',10,'Color','m');
                    text(tPts(stripeJumps(2)),1.25*pi,'jump #2','FontSize',10,'Color','m');
                end
                if hc == 2 & trialID == length(cond{condID}.allFlyData{flyID}.(trialType))
                    xlabel('time (s)');
                else
                    set(gca,'XTick',[]);
                end
                ylabel('position (rad)','Color',TCol);
                set(gca,'YTick',linspace(-pi,pi,5),'YTickLabels',{'-\pi','-\pi/2', '0', '\pi/2', '\pi'});
                colorbar;
                
                %% Plot OL
                figure(actPlotOL);
                % Plot the dark periods
                subplot(5,1,trialID+2*(hc-1));
                hold on;
                imagesc(tPts(OLPer),linspace(-pi,pi,size(DF,1)),DF(:,OLPer));
                plot(tPts(OLPer),stripePosPlt(OLPer),'b');
                plot(tPts(OLPer),headingPlt(OLPer),'k');
                axis tight;
                caxis([0 1]);
                if hc == 1 & trialID == 1
                    text(tPts(OLPer(1)),1.5*pi,strcat(cond{condID}.name,' - fly #',num2str(flyID)),'FontSize',14);
                end
                if hc == 2 & trialID == length(cond{condID}.allFlyData{flyID}.(trialType))
                    xlabel('time (s)');
                else
                    set(gca,'XTick',[]);
                end
                ylabel('position (rad)','Color',TCol);
                set(gca,'YTick',linspace(-pi,pi,5),'YTickLabels',{'-\pi','-\pi/2', '0', '\pi/2', '\pi'});
                colorbar;
                
                
                %% Plot direction CL
                figure(actPlotDirCL);
                % Plot the dark periods
                subplot(5,1,trialID+2*(hc-1));
                hold on;
                plot(tPts(darkPer),headingPlt(darkPer),'k');
                plot(tPts(CLPer),stripePosPlt(CLPer),'b');
                plot(tPts(tPtsDir),dirsPlt(plotPerDir),'g');
                
                xlim( [ tPts(1) max(tPts(plotPer)) ] )
                ylim([-pi pi])
                for jump = 1:length(stripeJumps);
                    line([tPts(stripeJumps(jump)) tPts(stripeJumps(jump))],...
                        [-pi pi],'Color','m');
                end
                axis tight;
                caxis([0 1]);
                if hc == 1 & trialID == 1
                    text(0,1.5*pi,strcat(cond{condID}.name,' - fly #',num2str(flyID)),'FontSize',14);
                    text(tPts(stripeJumps(1)),1.25*pi,'jump #1','FontSize',10,'Color','m');
                    text(tPts(stripeJumps(2)),1.25*pi,'jump #2','FontSize',10,'Color','m');
                end
                if hc == 2 & trialID == length(cond{condID}.allFlyData{flyID}.(trialType))
                    xlabel('time (s)');
                else
                    set(gca,'XTick',[]);
                end
                ylabel('position (rad)','Color',TCol);
                set(gca,'YTick',linspace(-pi,pi,5),'YTickLabels',{'-\pi','-\pi/2', '0', '\pi/2', '\pi'});    
                
                %% Plot vRot
                figure(actPlotvRot);
                % Plot the dark periods
                subplot(5,1,trialID+2*(hc-1));
                hold on;
                plot(tPts(darkPer),vRot(darkPer),'k');
                plot(tPts(CLPer),vRot(CLPer),'b');
                plot(tPts(tPtsDir),vDF(plotPerDir),'g');
                
                ymin = min( [min(vRot(darkPer)), min(vRot(CLPer))] );
                ymax = max( [max(vRot(darkPer)), max(vRot(CLPer))] );
                xlim( [ tPts(1) max(tPts(plotPer)) ] )

                for jump = 1:length(stripeJumps);
                    line([tPts(stripeJumps(jump)) tPts(stripeJumps(jump))],...
                        [ymin ymax],'Color','m');
                end
                axis tight;
                caxis([0 1]);
                if hc == 1 & trialID == 1
                    text(0,1.5*pi,strcat(cond{condID}.name,' - fly #',num2str(flyID)),'FontSize',14);
                    text(tPts(stripeJumps(1)),1.25*pi,'jump #1','FontSize',10,'Color','m');
                    text(tPts(stripeJumps(2)),1.25*pi,'jump #2','FontSize',10,'Color','m');
                end
                if hc == 2 & trialID == length(cond{condID}.allFlyData{flyID}.(trialType))
                    xlabel('time (s)');
                else
                    set(gca,'XTick',[]);
                end
                ylabel('velocity (rad/s)','Color',TCol);
                
                ylim( [ymin ymax] )
                %set(gca,'YTick',linspace(-pi,pi,5),'YTickLabels',{'-\pi','-\pi/2', '0', '\pi/2', '\pi'});
                
            end
        end
        
        figure(actPlotCL);
        colormap(brewermap(64, 'Blues'));
        set(actPlotCL,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
        print(actPlotCL,strcat(dir, 'raw_plots/',...
            cond{condID}.name(~isspace(cond{condID}.name)),'_Fly',num2str(flyID),'_ClosedLoop'),'-dpdf');
        delete(actPlotCL);

        figure(actPlotOL);
        colormap(brewermap(64, 'Blues'));
        set(actPlotOL,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
        print(actPlotOL,strcat(dir, 'raw_plots/',...
            cond{condID}.name(~isspace(cond{condID}.name)),'_Fly',num2str(flyID),'_OpenLoop'),'-dpdf');
        delete(actPlotOL);
        
        figure(actPlotDirCL);
        colormap(brewermap(64, 'Blues'));
        set(actPlotDirCL,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
        print(actPlotDirCL,strcat(dir, 'raw_plots/',...
            cond{condID}.name(~isspace(cond{condID}.name)),'_Fly',num2str(flyID),'_DirCL'),'-dpdf');
        delete(actPlotDirCL);
        
        figure(actPlotvRot);
        colormap(brewermap(64, 'Blues'));
        set(actPlotvRot,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
        print(actPlotvRot,strcat(dir, 'raw_plots/',...
            cond{condID}.name(~isspace(cond{condID}.name)),'_Fly',num2str(flyID),'_vRot'),'-dpdf');
        delete(actPlotvRot);
    end
end