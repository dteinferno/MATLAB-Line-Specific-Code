
function [tPts,darkPer,OLPer,CLPer,CWPer,CCWPer,DF,heading,headingPlt,vRot,vF,stripePos,stripePosPlt,stripeJumps]...
    = extractShiData(cond, flyID, trialType, trialID, smoothF, smoothV)
%extracts some of the more useful data from a container with a shibire
%experiment dataset. This is basically copied from dans 'LookAtShi'
    %cond is container
    %flyID is flynumber
    %trialtype is 'All' or 'All_30C'
    %trialID is 1-2 for all, 3-5 for 30C (although only use 4-5)
    %smooth specifies whether to apply an SG filter to the fluorescence
    %data (smoothF) and/or velocity data (smoothV)
    
    condID = 1; %I like separating my analysis by experiment, so only have condID == 1

    % Sort out the closed loop and dark periods
    [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch);
    darkEnd = find(diff(darkPer)>1);
    darkPer(darkEnd:end) = [];

    % Pull out the relevant data
    tPts = cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.OffsetRotMatch(:,1);
    tPts = tPts - tPts(1);
    heading = pi/180*cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.PosRotMatch;
    stripePos = cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.OffsetRotMatch(:,2);
    
    try
        DF = cond{condID}.allFlyData{flyID}.(trialType){trialID}.ROIaveMax-1;
    catch
        DF = cond{condID}.allFlyData{flyID}.(trialType){trialID}.GROIaveMax-1;
    end
    
    vRot = cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.vRot;
    vF = cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.vF;
    
    sgolayOrder = 3; %fit third order polynomial; could make this an argument as well 
    if smoothF > 0 %smoothen fluorescence
        sgolayFrames = smoothF;
        s = size(DF);
        new = [];
        for ind = 1:s(1);
            g = DF(ind,:);
            new(ind,:) = sgolayfilt(g,sgolayOrder,sgolayFrames);
        end
        DF = new;
    end
    if smoothV > 0 %smoothen velocity
        sgolayFrames = smoothV;
        vRot = sgolayfilt(vRot,sgolayOrder,sgolayFrames);
        vF = sgolayfilt(vF,sgolayOrder,sgolayFrames);
        heading = sgolayfilt(heading,sgolayOrder,sgolayFrames);
        stripePos = sgolayfilt(stripePos,sgolayOrder,sgolayFrames);
    end
    
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
    stripeJumps = find(diff(cond{condID}.allFlyData{flyID}.(trialType){trialID}.positionDatMatch.Trans) ~= 0);
    try
        stripeJumps(end) = [];
    catch
    end
                
                
                