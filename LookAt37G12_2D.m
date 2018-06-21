%% Clear out old data
clear;
clc;

cond = FlyDatLoad(1);

%% Plot the heading vs. the stripe position
% Find the PVA

condID = 1;
flyID = 1;

for trialID = 1:length(cond{condID}.allFlyData{flyID}.Stripe)

    FBAct = cond{condID}.allFlyData{flyID}.Stripe{trialID}.ROIaveMax-1;
    
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
    tAll = cond{condID}.allFlyData{flyID}.Stripe{trialID}.positionDatMatch.OffsetRotMatch(:,1);
    stripePos = cond{condID}.allFlyData{flyID}.Stripe{trialID}.positionDatMatch.OffsetRotMatch(:,2);
    
    subplot(3,1,trialID)
    hold on;
    plot(tAll,stripePos,'b');
    plot(tAll,meanAngRaw,'g');
end

%% Plot the heading vs. the world or cylinder position
% Find the PVA

condID = 1;
flyID = 1;

for trialID = 1:length(cond{condID}.allFlyData{flyID}.Cyl)

    FBAct = cond{condID}.allFlyData{flyID}.Cyl{trialID}.ROIaveMax-1;
    
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

    % Find the world position
    tAll = cond{condID}.allFlyData{flyID}.Cyl{trialID}.positionDatMatch.OffsetRotMatch(:,1);
    worldPos = cond{condID}.allFlyData{flyID}.Cyl{trialID}.positionDatMatch.OffsetRotMatch(:,2);
    
    % Find the cylinder position
    forPos = cond{condID}.allFlyData{flyID}.Cyl{trialID}.positionDatMatch.OffsetForMatch;
    latPos = cond{condID}.allFlyData{flyID}.Cyl{trialID}.positionDatMatch.OffsetLatMatch;
    cylAng = atan2(forPos,latPos)+worldPos+pi/2;
    
    subplot(3,2,trialID)
    hold on;
    plot(tAll,worldPos,'b');
    plot(tAll,cylAng,'k');
    plot(tAll,meanAngRaw,'g');
end