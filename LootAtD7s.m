%% Clear out old data
clear;
clc;
cond = FlyDatLoad(2);

%% Align the peaks and plot the normalized mean values as sorted by velocity

glomShift = 3;
vRMin = 0;
vRMax = 720;
vRSpan = 30;

redSpan = [2:17];
greenSpan = [2:17];

blues = brewermap(64, 'Blues');
greens(:,1) = blues(:,1);
greens(:,2) = blues(:,3);
greens(:,3) = blues(:,2);
magentas = blues;
magentas(:,1) = blues(:,3);
magentas(:,2) = blues(:,1);
magentas(:,3) = blues(:,3);
colorAssign = linspace(0,0.5,64);


[LPB1.dark.dataR, RPB1.dark.dataR, LPB1.dark.dataG, RPB1.dark.dataG] = ...
    BumpAlignOnePeak(cond{1}.allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'Dark',glomShift,0);
[LPB2.dark.dataR, RPB2.dark.dataR, LPB2.dark.dataG, RPB2.dark.dataG] = ...
    BumpAlignOnePeak(cond{2}.allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'Dark',glomShift,1);

vRBinNum = round((vRMax-vRMin)/vRSpan);

for condNow = 1:2
    
    if condNow == 1
        LPB = LPB1;
        RPB = RPB1;
    elseif condNow == 2
        LPB = LPB2;
        RPB = RPB2;
    end
    
    PBFig = figure('units','normalized','outerposition',[0 0 1 1]);

    num_ROIs = 8;
    angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
    angsraw = angsraw';

    plotRange = 180/vRSpan;
    for flyID = 1:cond{condNow}.numFlies

        LPBDataRAllStop = LPB.dark.dataR{flyID}.Stop;
        RPBDataRAllStop = RPB.dark.dataR{flyID}.Stop;
        LPBDataGAllStop = LPB.dark.dataG{flyID}.Stop;
        RPBDataGAllStop = RPB.dark.dataG{flyID}.Stop;

        % Plot the dark data
        LRStop = mean(LPBDataRAllStop,2);
        LRStop = (LRStop-min(LRStop(redSpan(1:8))))./(max(LRStop)-min(LRStop(redSpan(1:8))));
        RRStop = mean(RPBDataRAllStop,2);
        RRStop = (RRStop-min(RRStop(greenSpan(1:8))))./(max(RRStop)-min(RRStop(greenSpan(1:8))));
        RStop = vertcat(LRStop,RRStop);

        LGStop = mean(LPBDataGAllStop,2);
        LGStop = (LGStop-min(LGStop(greenSpan(1:8))))./(max(LGStop)-min(LGStop(greenSpan(1:8))));
        RGStop = median(RPBDataGAllStop,2);
        RGStop = (RGStop-min(RGStop(redSpan(1:8))))./(max(RGStop)-min(RGStop(redSpan(1:8))));
        GStop = vertcat(LGStop,RGStop);

        LRPk = find(LRStop == max(LRStop));
        LGPk = find(LGStop == max(LGStop));
        pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
        RRPk = find(RRStop == max(RRStop));
        RGPk = find(RGStop == max(RGStop));
        pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));

        subplot(2*cond{condNow}.numFlies,2*(plotRange+1),1+4*(plotRange+1)*(flyID-1))
        LBinR = mean(LPBDataRAllStop,2);
        RBinR = mean(RPBDataRAllStop,2);
        RBinR = circshift(RBinR,-1);
        actVals = vertcat(LBinR,RBinR)'-1;
        colorIm = zeros(2,18,3);
        for colorStep = 1:18
            colorNow = find(actVals(colorStep)<colorAssign)';
            colorNow(end+1) = 64;
            colorIm(1,colorStep,:) = magentas(colorNow(1),:);
        end
        LBinG =  mean(LPBDataGAllStop,2);
        RBinG = mean(RPBDataGAllStop,2);
        RBinG = circshift(RBinG,-1);
        actVals = vertcat(LBinG,RBinG)'-1;
        for colorStep = 1:18
            colorNow = find(actVals(colorStep)<colorAssign)';
            colorNow(end+1) = 64;
            colorIm(2,colorStep,:) = greens(colorNow(1),:);
        end
        image(colorIm);
        line([circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
            circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
            [0.5 1.5],'Color','k','LineWidth',2);
        line([circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
            circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
            [1.5 2.5],'Color','k','LineWidth',2);
%         text(2,0,num2str(pkDiffL));
%         text(11,0,num2str(pkDiffR));
        axis off;
        title(strcat('stopped, Fly ',num2str(flyID)));
        
        subplot(2*cond{condNow}.numFlies,2*(plotRange+1),2+4*(plotRange+1)*(flyID-1))
        RBinR = mean(RPBDataRAllStop,2);
        actVals = RBinR(2:9)-1;
        colorIm = zeros(2,8,3);
        for colorStep = 1:8
            colorNow = find(actVals(colorStep)<colorAssign)';
            colorNow(end+1) = 64;
            colorIm(1,colorStep,:) = magentas(colorNow(1),:);
        end
        RBinG = mean(RPBDataGAllStop,2);
        actVals = RBinG(2:9)-1;
        for colorStep = 1:8
            colorNow = find(actVals(colorStep)<colorAssign)';
            colorNow(end+1) = 64;
            colorIm(2,colorStep,:) = greens(colorNow(1),:);
        end
        image(colorIm);
        line([circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
            circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
            [0.5 1.5],'Color','k','LineWidth',2);
        line([circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
            circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
            [1.5 2.5],'Color','k','LineWidth',2);
%         text(2,0,num2str(pkDiffL));
%         text(11,0,num2str(pkDiffR));
        axis off;
        title(strcat('stopped, Fly ',num2str(flyID)));


        for binID=1:plotRange
            LPBDataRAllCW = LPB.dark.dataR{flyID}.CW{binID};
            RPBDataRAllCW = RPB.dark.dataR{flyID}.CW{binID};
            LPBDataGAllCW = LPB.dark.dataG{flyID}.CW{binID};
            RPBDataGAllCW = RPB.dark.dataG{flyID}.CW{binID};

            if ~isempty(LPBDataRAllCW)
                LRCW = mean(LPBDataRAllCW,2);
                LRCW = (LRCW-min(LRCW(redSpan(1:8))))./(max(LRCW)-min(LRCW(redSpan(1:8))));
                RRCW = mean(RPBDataRAllCW,2);
                RRCW = (RRCW-min(RRCW(greenSpan(1:8))))./(max(RRCW)-min(RRCW(greenSpan(1:8))));
                RCW = vertcat(LRCW,RRCW);
                LGCW = mean(LPBDataGAllCW,2);
                LGCW = (LGCW-min(LGCW(greenSpan(1:8))))./(max(LGCW)-min(LGCW(greenSpan(1:8))));
                RGCW = mean(RPBDataGAllCW,2);
                RGCW = (RGCW-min(RGCW(redSpan(1:8))))./(max(RGCW)-min(RGCW(redSpan(1:8))));
                GCW = vertcat(LGCW,RGCW);

                LRPk = find(LRCW == max(LRCW));
                LGPk = find(LGCW == max(LGCW));
                pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
                RRPk = find(RRCW == max(RRCW));
                RGPk = find(RGCW == max(RGCW));
                pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));


                subplot(2*cond{condNow}.numFlies,2*(plotRange+1),1+2*binID+4*(plotRange+1)*(flyID-1))
                LBinR = mean(LPBDataRAllCW,2);
                actVals = LBinR(2:9)-1;
                colorIm = zeros(2,8,3);
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(1,colorStep,:) = magentas(colorNow(1),:);
                end
                LBinG = mean(LPBDataGAllCW,2);
                actVals = LBinG(2:9)-1;
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(2,colorStep,:) = greens(colorNow(1),:);
                end
                image(colorIm);
                line([circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [0.5 1.5],'Color','k','LineWidth',2);
                line([circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [1.5 2.5],'Color','k','LineWidth',2);
%                 text(2,0,num2str(pkDiffL));
%                 text(11,0,num2str(pkDiffR));
                axis off;
                title(strcat(num2str(vRSpan*(binID-1)),'-',num2str(vRSpan*binID),' deg/s CW'));
                
                subplot(2*cond{condNow}.numFlies,2*(plotRange+1),2+2*binID+4*(plotRange+1)*(flyID-1))
                RBinR = mean(RPBDataRAllCW,2);
                actVals = RBinR(2:9)-1;
                colorIm = zeros(2,8,3);
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(1,colorStep,:) = magentas(colorNow(1),:);
                end
                RBinG = mean(RPBDataGAllCW,2);
                actVals = RBinG(2:9)-1;
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(2,colorStep,:) = greens(colorNow(1),:);
                end
                image(colorIm);
                line([circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [0.5 1.5],'Color','k','LineWidth',2);
                line([circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [1.5 2.5],'Color','k','LineWidth',2);
%                 text(2,0,num2str(pkDiffL));
%                 text(11,0,num2str(pkDiffR));
                axis off;
                title(strcat(num2str(vRSpan*(binID-1)),'-',num2str(vRSpan*binID),' deg/s CW'));
            end

            LPBDataRAllCCW = LPB.dark.dataR{flyID}.CCW{binID};
            RPBDataRAllCCW = RPB.dark.dataR{flyID}.CCW{binID};
            LPBDataGAllCCW = LPB.dark.dataG{flyID}.CCW{binID};
            RPBDataGAllCCW = RPB.dark.dataG{flyID}.CCW{binID};

            if ~isempty(LPBDataRAllCCW)
                LRCCW = mean(LPBDataRAllCCW,2);
                LRCCW = (LRCCW-min(LRCCW(redSpan(1:8))))./(max(LRCCW)-min(LRCCW(redSpan(1:8))));
                RRCCW = mean(RPBDataRAllCCW,2);
                RRCCW = (RRCCW-min(RRCCW(greenSpan(1:8))))./(max(RRCCW)-min(RRCCW(greenSpan(1:8))));
                RCCW = vertcat(LRCCW,RRCCW);
                LGCCW = mean(LPBDataGAllCCW,2);
                LGCCW = (LGCCW-min(LGCCW(greenSpan(1:8))))./(max(LGCCW)-min(LGCCW(greenSpan(1:8))));
                RGCCW = mean(RPBDataGAllCCW,2);
                RGCCW = (RGCCW-min(RGCCW(redSpan(1:8))))./(max(RGCCW)-min(RGCCW(redSpan(1:8))));
                GCCW = vertcat(LGCCW,RGCCW);

                LRPk = find(LRCCW == max(LRCCW));
                LGPk = find(LGCCW == max(LGCCW));
                pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
                RRPk = find(RRCCW == max(RRCCW));
                RGPk = find(RGCCW == max(RGCCW));
                pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));

                subplot(2*cond{condNow}.numFlies,2*(plotRange+1),1+2*binID+4*(plotRange+1)*(flyID-1)+2*(plotRange+1))
                LBinR = mean(LPBDataRAllCCW,2);
                actVals = LBinR(2:9)-1;
                colorIm = zeros(2,8,3);
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(1,colorStep,:) = magentas(colorNow(1),:);
                end
                LBinG = mean(LPBDataGAllCCW,2);
                actVals = LBinG(2:9)-1;
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(2,colorStep,:) = greens(colorNow(1),:);
                end
                image(colorIm);
                line([circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [0.5 1.5],'Color','k','LineWidth',2);
                line([circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [1.5 2.5],'Color','k','LineWidth',2);
%                 text(2,0,num2str(pkDiffL));
%                 text(11,0,num2str(pkDiffR));
                axis off;
                title(strcat(num2str(vRSpan*(binID-1)),'-',num2str(vRSpan*binID),' deg/s CCW'));
                
                subplot(2*cond{condNow}.numFlies,2*(plotRange+1),2+2*binID+4*(plotRange+1)*(flyID-1)+2*(plotRange+1))
                RBinR = mean(RPBDataRAllCCW,2);
                actVals = RBinR(2:9)-1;
                colorIm = zeros(2,8,3);
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(1,colorStep,:) = magentas(colorNow(1),:);
                end
                RBinG = mean(RPBDataGAllCCW,2);
                actVals = RBinG(2:9)-1;
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(2,colorStep,:) = greens(colorNow(1),:);
                end
                image(colorIm);
                line([circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [0.5 1.5],'Color','k','LineWidth',2);
                line([circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [1.5 2.5],'Color','k','LineWidth',2);
%                 text(2,0,num2str(pkDiffL));
%                 text(11,0,num2str(pkDiffR));
                axis off;
                title(strcat(num2str(vRSpan*(binID-1)),'-',num2str(vRSpan*binID),' deg/s CCW'));
            end
        end
    end

%     figure(PBFig);
    [ax,h] = suplabel(cond{condNow}.name,'t',[0.1 0.1 0.85 0.85]);
    % set(PBFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    % print(PBFig,strcat('All_vRotVsPBPosAll_OnePeak_NormColor_RAlign_greens'),'-dpdf');
end

%% Plot the offsets as a function of velocity

glomShift = 3;
vRMin = 0;
vRMax = 720;
vRSpan = 30;

redSpan = [2:17];
greenSpan = [2:17];

blues = brewermap(64, 'Blues');
greens(:,1) = blues(:,1);
greens(:,2) = blues(:,3);
greens(:,3) = blues(:,2);
magentas = blues;
magentas(:,1) = blues(:,3);
magentas(:,2) = blues(:,1);
magentas(:,3) = blues(:,3);
colorAssign = linspace(0,0.5,64);


[LPB1.dark.dataR, RPB1.dark.dataR, LPB1.dark.dataG, RPB1.dark.dataG] = ...
    BumpAlignOnePeak(cond{1}.allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'Dark',glomShift,0);
[LPB2.dark.dataR, RPB2.dark.dataR, LPB2.dark.dataG, RPB2.dark.dataG] = ...
    BumpAlignOnePeak(cond{2}.allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'Dark',glomShift,1);

vRBinNum = round((vRMax-vRMin)/vRSpan);

for condNow = 1:2
    
    if condNow == 1
        LPB = LPB1;
        RPB = RPB1;
    elseif condNow == 2
        LPB = LPB2;
        RPB = RPB2;
    end
    
    PBFig = figure('units','normalized','outerposition',[0 0 1 1]);

    num_ROIs = 8;
    angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
    angsraw = angsraw';

    plotRange = 180/vRSpan;
    for flyID = 1:cond{condNow}.numFlies

        LPBDataRAllStop = LPB.dark.dataR{flyID}.Stop;
        RPBDataRAllStop = RPB.dark.dataR{flyID}.Stop;
        LPBDataGAllStop = LPB.dark.dataG{flyID}.Stop;
        RPBDataGAllStop = RPB.dark.dataG{flyID}.Stop;

        % Plot the dark data
        LRStop = mean(LPBDataRAllStop,2);
        LRStop = (LRStop-min(LRStop(redSpan(1:8))))./(max(LRStop)-min(LRStop(redSpan(1:8))));
        RRStop = mean(RPBDataRAllStop,2);
        RRStop = (RRStop-min(RRStop(greenSpan(1:8))))./(max(RRStop)-min(RRStop(greenSpan(1:8))));
        RStop = vertcat(LRStop,RRStop);

        LGStop = mean(LPBDataGAllStop,2);
        LGStop = (LGStop-min(LGStop(greenSpan(1:8))))./(max(LGStop)-min(LGStop(greenSpan(1:8))));
        RGStop = median(RPBDataGAllStop,2);
        RGStop = (RGStop-min(RGStop(redSpan(1:8))))./(max(RGStop)-min(RGStop(redSpan(1:8))));
        GStop = vertcat(LGStop,RGStop);

        LRPk = find(LRStop == max(LRStop));
        LGPk = find(LGStop == max(LGStop));
        pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
        RRPk = find(RRStop == max(RRStop));
        RGPk = find(RGStop == max(RGStop));
        pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));
        
        PVADiffL = abs(circ_mean(angsraw,LGStop(2:9))-circ_mean(angsraw,LRStop(2:9)))...
            *num_ROIs/(2*pi);
        PVADiffL = min(PVADiffL,8-PVADiffL);
        PVADiffR = abs(circ_mean(angsraw,RGStop(2:9))-circ_mean(angsraw,RRStop(2:9)))...
            *num_ROIs/(2*pi);
        PVADiffR = min(PVADiffR,8-PVADiffR);

        subplot(2,2,1);
        hold on;
        scatter(0,pkDiffL);
        
        subplot(2,2,2);
        hold on;
        scatter(0,pkDiffR);
        
        subplot(2,2,3);
        hold on;
        scatter(0,PVADiffL);
        
        subplot(2,2,4);
        hold on;
        scatter(0,PVADiffR);

        for binID=1:plotRange
            LPBDataRAllCW = LPB.dark.dataR{flyID}.CW{binID};
            RPBDataRAllCW = RPB.dark.dataR{flyID}.CW{binID};
            LPBDataGAllCW = LPB.dark.dataG{flyID}.CW{binID};
            RPBDataGAllCW = RPB.dark.dataG{flyID}.CW{binID};

            if ~isempty(LPBDataRAllCW)
                LRCW = mean(LPBDataRAllCW,2);
                LRCW = (LRCW-min(LRCW(redSpan(1:8))))./(max(LRCW)-min(LRCW(redSpan(1:8))));
                RRCW = mean(RPBDataRAllCW,2);
                RRCW = (RRCW-min(RRCW(greenSpan(1:8))))./(max(RRCW)-min(RRCW(greenSpan(1:8))));
                RCW = vertcat(LRCW,RRCW);
                LGCW = mean(LPBDataGAllCW,2);
                LGCW = (LGCW-min(LGCW(greenSpan(1:8))))./(max(LGCW)-min(LGCW(greenSpan(1:8))));
                RGCW = mean(RPBDataGAllCW,2);
                RGCW = (RGCW-min(RGCW(redSpan(1:8))))./(max(RGCW)-min(RGCW(redSpan(1:8))));
                GCW = vertcat(LGCW,RGCW);

                LRPk = find(LRCW == max(LRCW));
                LGPk = find(LGCW == max(LGCW));
                pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
                RRPk = find(RRCW == max(RRCW));
                RGPk = find(RGCW == max(RGCW));
                pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));
                
                PVADiffL = abs(circ_mean(angsraw,LGCW(2:9))-circ_mean(angsraw,LRCW(2:9)))...
                    *num_ROIs/(2*pi);
                PVADiffL = min(PVADiffL,8-PVADiffL);
                PVADiffR = abs(circ_mean(angsraw,RGCW(2:9))-circ_mean(angsraw,RRCW(2:9)))...
                    *num_ROIs/(2*pi);
                PVADiffR = min(PVADiffR,8-PVADiffR);
                
                subplot(2,2,1);
                scatter(binID,pkDiffL);
                
                subplot(2,2,2);
                scatter(binID,pkDiffR);
                
                subplot(2,2,3);
                scatter(binID,PVADiffL);
                
                subplot(2,2,4);
                scatter(binID,PVADiffR);

            end

            LPBDataRAllCCW = LPB.dark.dataR{flyID}.CCW{binID};
            RPBDataRAllCCW = RPB.dark.dataR{flyID}.CCW{binID};
            LPBDataGAllCCW = LPB.dark.dataG{flyID}.CCW{binID};
            RPBDataGAllCCW = RPB.dark.dataG{flyID}.CCW{binID};

            if ~isempty(LPBDataRAllCCW)
                LRCCW = mean(LPBDataRAllCCW,2);
                LRCCW = (LRCCW-min(LRCCW(redSpan(1:8))))./(max(LRCCW)-min(LRCCW(redSpan(1:8))));
                RRCCW = mean(RPBDataRAllCCW,2);
                RRCCW = (RRCCW-min(RRCCW(greenSpan(1:8))))./(max(RRCCW)-min(RRCCW(greenSpan(1:8))));
                RCCW = vertcat(LRCCW,RRCCW);
                LGCCW = mean(LPBDataGAllCCW,2);
                LGCCW = (LGCCW-min(LGCCW(greenSpan(1:8))))./(max(LGCCW)-min(LGCCW(greenSpan(1:8))));
                RGCCW = mean(RPBDataGAllCCW,2);
                RGCCW = (RGCCW-min(RGCCW(redSpan(1:8))))./(max(RGCCW)-min(RGCCW(redSpan(1:8))));
                GCCW = vertcat(LGCCW,RGCCW);

                LRPk = find(LRCCW == max(LRCCW));
                LGPk = find(LGCCW == max(LGCCW));
                pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
                RRPk = find(RRCCW == max(RRCCW));
                RGPk = find(RGCCW == max(RGCCW));
                pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));
                
                PVADiffL = abs(circ_mean(angsraw,LGCCW(2:9))-circ_mean(angsraw,LRCCW(2:9)))...
                    *num_ROIs/(2*pi);
                PVADiffL = min(PVADiffL,8-PVADiffL);
                PVADiffR = abs(circ_mean(angsraw,RGCCW(2:9))-circ_mean(angsraw,RRCCW(2:9)))...
                    *num_ROIs/(2*pi);
                PVADiffR = min(PVADiffR,8-PVADiffR);

                subplot(2,2,1);
                scatter(-binID,pkDiffL);
                
                subplot(2,2,2);
                scatter(-binID,pkDiffR);
                
                subplot(2,2,3);
                scatter(-binID,PVADiffL);
                
                subplot(2,2,4);
                scatter(-binID,PVADiffR);
            end
        end
    end

%     figure(PBFig);
    [ax,h] = suplabel(cond{condNow}.name,'t',[0.1 0.1 0.85 0.85]);
    % set(PBFig,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
    % print(PBFig,strcat('All_vRotVsPBPosAll_OnePeak_NormColor_RAlign_greens'),'-dpdf');
end

%% Do LM with activity peak, rotational velocity, and forward velocity

PBRange = [2:9];

numtPts = 51;

for condNow = 1:2
    
    figure;
    
    for flyID = 1:cond{condNow}.numFlies
        
        lrCoefR = zeros(numtPts,3);
        lrCoefG = zeros(numtPts,3);
        
        vRot = [];
        vF = [];
        RMaxVals = [];
        GMaxVals = [];
        
        for trial = 1:length(cond{condNow}.allFlyData{flyID}.Dark)
           vRot = vertcat(vRot, ...
               cond{condNow}.allFlyData{flyID}.Dark{trial}.positionDatMatch.vRot);
           vF = vertcat(vF, ...
               cond{condNow}.allFlyData{flyID}.Dark{trial}.positionDatMatch.vF);
           RMaxVals = vertcat(RMaxVals,...
               max(cond{condNow}.allFlyData{flyID}.Dark{trial}.RROIaveMax(PBRange,2:end))');
           GMaxVals = vertcat(GMaxVals,...
               max(cond{condNow}.allFlyData{flyID}.Dark{trial}.GROIaveMax(PBRange,2:end))');
        end
        
        vCW = zeros(size(vRot));
        vCCW = zeros(size(vRot));
        vPos = find(vRot > 0);
        vNeg = find(vRot < 0);
        vCW(vPos) = vRot(vPos);
        vCCW(vNeg) = abs(vRot(vNeg));
        vPred = horzcat(zscore(vCW),zscore(vCCW),zscore(vF)); 
        
        for tStep = 1:numtPts
            lrCoefR(tStep,:) = ...
                regress(...
                zscore(RMaxVals(tStep:end-numtPts+tStep)),...
                vPred(floor(numtPts/2):end-ceil(numtPts/2),:));
            lrCoefG(tStep,:) = ...
                regress(...
                zscore(GMaxVals(tStep:end-numtPts+tStep)),...
                vPred(floor(numtPts/2):end-ceil(numtPts/2),:));
        end
        
        subplot(2,cond{condNow}.numFlies,flyID);
        hold on;
        plot(lrCoefR(:,1),'b');
        plot(lrCoefR(:,2),'r');
        plot(lrCoefR(:,3),'m');
        
        subplot(2,cond{condNow}.numFlies,flyID+cond{condNow}.numFlies);
        hold on;
        plot(lrCoefG(:,1),'b');
        plot(lrCoefG(:,2),'r');
        plot(lrCoefG(:,3),'m');
        
    end
    
end

%% Look at relationship between red and green peak values 

PBRange = [2:9];

for condNow = 1:2
    
    figure;
    
    for flyID = 1:cond{condNow}.numFlies
        
        RMaxVals = [];
        GMaxVals = [];
        
        for trial = 1:length(cond{condNow}.allFlyData{flyID}.Dark)
           RMaxVals = vertcat(RMaxVals,...
               max(cond{condNow}.allFlyData{flyID}.Dark{trial}.RROIaveMax(PBRange,:))');
           GMaxVals = vertcat(GMaxVals,...
               max(cond{condNow}.allFlyData{flyID}.Dark{trial}.GROIaveMax(PBRange,:))');
        end
        
        subplot(1,cond{condNow}.numFlies,flyID);
        scatter(RMaxVals,GMaxVals,20,'filled')
        alpha(0.1);
        
    end
    
end

%% Load an example tiff
pathName = 'C:\Users\turnerevansd\Documents\D7s\G60D05R55G08\20170717\';
tiffFilename = 'Fly2_3-4day_6fx60D05_jRGCx55G08_Dark_00001.tif';

behavFilename = 'Fly2_3-4day_6fx60D05_jRGCx55G08_Dark_01.TXT';
positionDat = VRDatLoad(behavFilename,pathName,0);

% Load the imaging stack and get the behavioral data
[RstackMaxInt, GstackMaxInt, RstackMean, GstackMean] = ImDatLoadBigtiff2Color(tiffFilename,pathName);

%% Preprocess the imaging stack
% Subtract the background
RMean = mean(RstackMaxInt,3);
GMean = mean(GstackMaxInt,3);
for frame = 1:size(RstackMaxInt,3)
    RstackSub(:,:,frame) = RstackMaxInt(:,:,frame) - RMean;
    GstackSub(:,:,frame) = GstackMaxInt(:,:,frame) - GMean;
end

% Gaussian filter the stacks
gaussianSize = [7 7];
gaussianSigma = 3;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);
RstackXYfiltMax = double(zeros(size(RstackMaxInt)));
GstackXYfiltMax = double(zeros(size(GstackMaxInt)));

h = waitbar(0.0,'Gaussian filtering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Gaussian filtering TIFF stack...');
for i = 1:size(RstackMaxInt,3)
    if mod(i,100)==0
        waitbar(i/length(RstackMaxInt),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackMaxInt))]);
    end
    RstackXYfiltMax(:,:,i) = imfilter(RstackMaxInt(:,:,i),Gxy,'replicate','conv');
    GstackXYfiltMax(:,:,i) = imfilter(GstackMaxInt(:,:,i),Gxy,'replicate','conv');
end
delete(h);

%% Make a movie showing the combined activity

num_planes = round(length(positionDat.tFrameGrab)/size(RstackMaxInt,3));
tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
maxFG = min(round(length(positionDat.tFrameGrab)/num_planes),size(RstackMaxInt,3));

% Find the min and max of the imaging stack and the imaging rate
RmaxCa = max(max(max(RstackXYfiltMax)));
GmaxCa = max(max(max(GstackXYfiltMax)));
RminCa = RmaxCa/4;% min(min(min(RstackXYfiltMax)));
GminCa = GmaxCa/7;% min(min(min(GstackXYfiltMax)));
imFreq = 10000/mean(diff(positionDat.tFrameGrab))/num_planes; 

fmStart = round(imFreq*35);
fmEnd = round(imFreq*40);

% Open the movie object and make the figure
writerObj = VideoWriter(strcat(pathName,'D7.avi'));
writerObj.FrameRate= 2*imFreq;
open(writerObj);
frameNow = figure('Position',[50 50 500 500],'Color','k');

% Generate the movie
for fm = fmStart:fmEnd

    % Create the colored images from the two stacks
    curRFrame = 2*(RstackXYfiltMax(60:160,:,fm+minFG-1)-...
        RminCa)./(RmaxCa-RminCa);
    curGFrame = 6*(GstackXYfiltMax(60:160,:,fm+minFG-1)-...
        GminCa)./(GmaxCa-GminCa);

    RIm = zeros([size(curRFrame) 3]);
    RIm(:,:,1) = curRFrame;
    RIm(:,:,3) = curRFrame;
    GIm = zeros([size(curGFrame) 3]);
    GIm(:,:,2) = curGFrame;

    subplot(3,1,1);
    cla;
    imshow(RIm);
    axis off;
    
    subplot(3,1,2);
    cla;
    imshow(GIm);
    axis off;
    
    subplot(3,1,3);
    cla;
    imshow(RIm+GIm);
    axis off;
  
    
    frame = getframe(frameNow);
    writeVideo(writerObj,frame);
end

delete(frameNow);
close(writerObj);

%% Show stills of the bump moving in both populations across three time points and the max. int. plot over time

% Find the min and max of the imaging stack and the imaging rate
RmaxCa = max(max(max(RstackXYfiltMax)));
GmaxCa = max(max(max(GstackXYfiltMax)));
RminCa = RmaxCa/4;% min(min(min(RstackXYfiltMax)));
GminCa = GmaxCa/6;% min(min(min(GstackXYfiltMax)));
imFreq = 10000/mean(diff(positionDat.tFrameGrab))/num_planes; 

exPoint1 = round(imFreq*38.7);
exPoint2 = round(imFreq*42.7);
exPoint3 = round(imFreq*51.7);

TwoPop = figure('Position',[50 50 800 1000]);

RMean = mean(RstackMaxInt(60:160,:,:),3);
RMean = RMean./max(max(RMean));

subplot(6,4,1);
imshow(RMean);

GMean = mean(GstackMaxInt(60:160,:,:),3);
GMean = GMean./max(max(GMean));

subplot(6,4,5);
imshow(GMean);

% Create the colored images from the two stacks
curRFrame = 2*(RstackXYfiltMax(60:160,:,exPoint1)-...
    RminCa)./(RmaxCa-RminCa);
curGFrame = 6*(GstackXYfiltMax(60:160,:,exPoint1)-...
    GminCa)./(GmaxCa-GminCa);

RIm = zeros([size(curRFrame) 3]);
RIm(:,:,1) = curRFrame;
RIm(:,:,3) = curRFrame;
GIm = zeros([size(curGFrame) 3]);
GIm(:,:,2) = curGFrame;

subplot(6,4,2);
cla;
imshow(RIm);
axis off;

subplot(6,4,6);
cla;
imshow(GIm);
axis off;

subplot(6,4,10);
cla;
imshow(GIm+RIm);
axis off;

colorBarIm = zeros(2,50,3);
colorBarIm(1,:,1) = linspace(1,50,50)./50;
colorBarIm(1,:,3) = linspace(1,50,50)./50;
colorBarIm(2,:,2) = linspace(1,50,50)./50;

subplot(6,4,9);
image(1:50,1:2,colorBarIm)
set(gca,'XTick',[1 50],'XTickLabels',...
    {num2str(GminCa) num2str((GmaxCa-GminCa)/6+GminCa)});
text(-5,0.25,num2str(RminCa));
text(45,0.25,num2str((RmaxCa-RminCa)/2+RminCa));

curRFrame = 2*(RstackXYfiltMax(60:160,:,exPoint2)-...
    RminCa)./(RmaxCa-RminCa);
curGFrame = 6*(GstackXYfiltMax(60:160,:,exPoint2)-...
    GminCa)./(GmaxCa-GminCa);

RIm = zeros([size(curRFrame) 3]);
RIm(:,:,1) = curRFrame;
RIm(:,:,3) = curRFrame;
GIm = zeros([size(curGFrame) 3]);
GIm(:,:,2) = curGFrame;

subplot(6,4,3);
cla;
imshow(RIm);
axis off;

subplot(6,4,7);
cla;
imshow(GIm);
axis off;

subplot(6,4,11);
cla;
imshow(GIm+RIm);
axis off;

curRFrame = 2*(RstackXYfiltMax(60:160,:,exPoint3)-...
    RminCa)./(RmaxCa-RminCa);
curGFrame = 6*(GstackXYfiltMax(60:160,:,exPoint3)-...
    GminCa)./(GmaxCa-GminCa);

RIm = zeros([size(curRFrame) 3]);
RIm(:,:,1) = curRFrame;
RIm(:,:,3) = curRFrame;
GIm = zeros([size(curGFrame) 3]);
GIm(:,:,2) = curGFrame;

subplot(6,4,4);
cla;
imshow(RIm);
axis off;

subplot(6,4,8);
cla;
imshow(GIm);
axis off;

subplot(6,4,12);
cla;
imshow(GIm+RIm);
axis off;

fmStart = round(imFreq*35);
fmEnd = round(imFreq*55);

tAll = cond{1}.allFlyData{2}.Dark{1}.positionDatMatch.OffsetRotMatch(fmStart:fmEnd,1);
t1 = cond{1}.allFlyData{2}.Dark{1}.positionDatMatch.OffsetRotMatch(exPoint1,1);
t2 = cond{1}.allFlyData{2}.Dark{1}.positionDatMatch.OffsetRotMatch(exPoint2,1);
t3 = cond{1}.allFlyData{2}.Dark{1}.positionDatMatch.OffsetRotMatch(exPoint3,1);
RData = cond{1}.allFlyData{2}.Dark{1}.RROIaveMax(2:17,fmStart:fmEnd);
GData = cond{1}.allFlyData{2}.Dark{1}.GROIaveMax(2:17,fmStart:fmEnd);

RIm = zeros([size(RData) 3]);
GIm = zeros([size(GData) 3]);

RMax = max(max(RData));
RMin = min(min(RData));
GMax = max(max(GData));
GMin = min(min(GData));

RIm(:,:,1) = (RData-RMin)./(RMax-RMin);
RIm(:,:,3) = (RData-RMin)./(RMax-RMin);
GIm(:,:,2) = (GData-GMin)./(GMax-GMin);

subplot(6,4,[14:16]);
image(tAll,2:17,RIm);
set(gca,'XTick',[]);
ylabel('PB glomerulus');
title('\Delta7 neurons');
line([t1 t1],[1 18],'Color','w','LineStyle','--');
line([t2 t2],[1 18],'Color','w','LineStyle','--');
line([t3 t3],[1 18],'Color','w','LineStyle','--');
line([tAll(1) tAll(end)],[9.5 9.5],'Color','w');
set(gca,'XColor','w');

subplot(6,4,[18:20]);
image(tAll,2:17,GIm);
ylabel('PB glomerulus');
title('compass neurons');
line([t1 t1],[1 18],'Color','w','LineStyle','--');
line([t2 t2],[1 18],'Color','w','LineStyle','--');
line([t3 t3],[1 18],'Color','w','LineStyle','--');
line([tAll(1) tAll(end)],[9.5 9.5],'Color','w');


subplot(6,4,[22:24]);
image(tAll,2:17,RIm+GIm);
xlabel('time (s)')
ylabel('PB glomerulus');
line([t1 t1],[1 18],'Color','w','LineStyle','--');
line([t2 t2],[1 18],'Color','w','LineStyle','--');
line([t3 t3],[1 18],'Color','w','LineStyle','--');
line([tAll(1) tAll(end)],[9.5 9.5],'Color','w');

set(TwoPop,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
print(TwoPop,strcat('C:\Users\turnerevansd\Documents\D7s\Results\ExampleBumps'),'-dpdf');

%% Look at relative bump widths of P-ENs and E-PG


