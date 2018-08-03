function [ominusR, oplusR, ominusL, oplusL, oAbsR, oAbsL, actgR, actrR, actgL, actrL, slopes] = plot_PB_Kris(dat, green, red, fly, trial, dir)
%given a single trial for a single fly in the PB, plots preliminary data: heading vs green/red directions, vR/vF vs.
%green and red intensity and vector sum magnitude (NOT normalized), vRot vs offset.
%Also makes boxplots of intensity & offset, and scatterplots of green/red offset vs vRot as
%well as the rate of change of the direction of the green/red vector sum vs. vRot.
%Does all of this for the 'complete' dataset and for thresholded data by
%the magnitude of the green and red vector sums.
%Does this for left and right PB

%% Extract and process data

vR = dat.positionDatMatch.vRot( dat.positionDatMatch.Closed(1:length(dat.positionDatMatch.vRot)) == 1 ) ; %CCW velocity
vF = dat.positionDatMatch.vF(dat.positionDatMatch.Closed(1:length(dat.positionDatMatch.vF))== 1); %forward velocity


ts = dat.positionDatMatch.OffsetRotMatch(:,1);
ts = ts(dat.positionDatMatch.Closed== 1);

L = length(dat.GROIaveMax);

heading = dat.positionDatMatch.PosRotMatch(1:L) * pi/180;
heading = heading(dat.positionDatMatch.Closed== 1);%CCW is positive

datG = flipud(dat.GROIaveMax)-1; %we're numbering from 1 at L9 to 18 at R9 on Tanya's diagram
datR = flipud(dat.RROIaveMax)-1; %this means CCW rotation (+ve vrot) gives higher numbers

s = size(datG);
s = s(1);

smooth = 3;
%should change to SG filtering
if smooth > 0
    display('smoothening data');
    ts = Smooth(ts, smooth);
    vR = Smooth(vR, smooth);
    heading = Smooth(heading, smooth);
    vF = Smooth(vF, smooth);
    newG = [];
    newG = [];
    for i = 1:s;
        g = datG(i,:);
        newG(i,:) = Smooth(g(dat.positionDatMatch.Closed== 1), smooth);
        r = datR(i,:);
        newR(i,:) = Smooth(r(dat.positionDatMatch.Closed== 1), smooth);
        %newR(i,:) = Smooth(datR(i,:), smooth);
    end
    datG = newG;
    datR = newR; 
end

[mGl, mGr, dGl, dGr] = getPBVec(datG); %get vector magnitudes and direction
[mRl, mRr, dRl, dRr] = getPBVec(datR);

mGt = mGr+mGl; %total intensities across PB
mRt = mRr+mRl;

L = length(vR); %number of data points we have velocity for

%%prune away stuff with too low intensity

prune = 0.20;

if prune > 0
    rlim = prune*mean(mRt);
    glim = prune*mean(mGt);

    tp = ts( mRt > rlim );
    vRp = vR( mRt(1:end-1) > rlim );
    dGrp = dGr( mRt > rlim );
    mGrp = mGr( mRt > rlim );
    dRrp = dRr( mRt > rlim );
    mRrp = mRr( mRt > rlim );
    dGlp = dGl( mRt > rlim );
    mGlp = mGl( mRt > rlim );
    dRlp = dRl( mRt > rlim );
    mRlp = mRl( mRt > rlim );
    
    mGpt = mGt( mRt > rlim );
    
    tp = tp( mGpt > glim );
    vRp = vRp( mGpt(1:end-1) > glim ); 
    dGrp = dGrp( mGpt > glim );
    dRrp = dRrp( mGpt > glim );
    mRrp = mRrp( mGpt > glim );
    mGrp = mGrp( mGpt > glim );
    dGlp = dGlp( mGpt > glim );
    dRlp = dRlp( mGpt > glim );
    mRlp = mRlp( mGpt > glim );
    mGlp = mGlp( mGpt > glim );
end


rGr = getRates(dGr, ts, 9);%rate of change of direction of Green vector sum. +ve is moving to the right in the PB (e.g. CCW fly rotation)
rRr = getRates(dRr, ts, 9);%rate of change of direction of Red vector sum
rGrp = getRates(dGrp, tp, 9);%pruned
rRrp = getRates(dRrp, tp, 9);
rGl = getRates(dGl, ts, 9);%left PB
rRl = getRates(dRl, ts, 9);
rGlp = getRates(dGlp, tp, 9);
rRlp = getRates(dRlp, tp, 9);

%cut off last value to amke compatible with our rate data
dGr = dGr(1:L);
dRr = dRr(1:L);
mGr = mGr(1:L);
mRr = mRr(1:L);
dGl = dGl(1:L);
dRl = dRl(1:L);
mGl = mGl(1:L);
mRl = mRl(1:L);
ts = ts(1:L);
heading = heading(1:L);

offsetr = getOffset(dGr, dRr, 9)'; %get offset of green relative to red
offsetrp = getOffset(dGrp, dRrp, 9)'; %pruned
offsetl = getOffset(dGl, dRl, 9)'; %left
offsetlp = getOffset(dGlp, dRlp, 9)';

%% plot raw data

% Plot the profiles
act = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

% Plot the rotational velocity
subplot(3,1,1);
hold on;

plot(ts,getCont(dGr, 5),'g'); %vectorPos
plot(ts,getCont(dRr, 5),'r');%vectorPoos
plot(ts,getCont(dGl, 5),'g'); %vectorPos

[hAx,hLine1,hLine2] = plotyy(ts,getCont(dRl, 5),ts,getCont(heading, pi));

set(hLine1,'Color','r');
set(hLine2,'Color','b');

plot([ts(1), ts(end)], [9, 9], 'k')
legend({green, red, green, red, 'G9', 'Heading'});
ylabel('position (glomerulus)');

ylim(hAx(2), [-pi 3*pi])
ylim(hAx(1), [0 18] );
xlim(hAx(1), [ ts(1) ts(end) ] )
xlim(hAx(2), [ ts(1) ts(end) ] )
set(hAx,{'ycolor'},{'k';'k'})


%%forward velocity
subplot(3,1,2);
absvF = abs(vF);
hold on
plot(ts,mGr,'g--'); %vectorMag
plot(ts,mGl,'g')
plot(ts,mRr,'r--')

[hAx,hLine1,hLine2] = plotyy(ts,mRl,ts,absvF);
set(hAx,{'ycolor'},{'k';'k'})

ylabel(hAx(1),'|Position Vector|') % left y-axis 
ylabel(hAx(2),'|Forward Velocity|') % right y-axis
xlim(hAx(1), [ ts(1) ts(end) ] )
xlim(hAx(2), [ ts(1) ts(end) ] )
ylim(hAx(1), [0 max( [ max(mRr) max(mGr) max(mRl) max(mGl) ] )])

ylim(hAx(2), [min( absvF ) max(absvF)])
set(hLine1,'Color','r');
set(hLine2,'Color','k');
legend({strcat(green, ' R'), strcat(green, ' L'), strcat(red, 'R'), strcat(red, 'L'), 'vF'});

%%rotational velocity
subplot(3,1,3);
absvR = abs(vR);
hold on
plot(ts,mGr,'g--'); %vectorMag
plot(ts,mGl,'g')
plot(ts,mRr,'r--')
[hAx,hLine1,hLine2] = plotyy(ts,mRl,ts,absvR);
set(hAx,{'ycolor'},{'k';'k'})

ylabel(hAx(1),'|Activity Vector|') % left y-axis 
ylabel(hAx(2),'|Position Velocity|') % right y-axis
xlim(hAx(1), [ ts(1) ts(end) ] )
xlim(hAx(2), [ ts(1) ts(end) ] )
ylim(hAx(1), [0 max( [ max(mRr) max(mGr) max(mRl) max(mGl) ] )])
ylim(hAx(2), [min( absvR ) max(absvR)])
set(hLine1,'Color','r');
set(hLine2,'Color','k');
legend({strcat(green, ' R'), strcat(green, ' L'), strcat(red, 'R'), strcat(red, 'L'), 'vRot'});

%% plot vs offset and velocity


act2 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

%rate of change
subplot(2,1,1);
ymax = max( max(vR), abs(min(vR)) );
ymax2 = max( [ max(abs(rGl)) max(abs(rRl)) max(abs(rGr)) max(abs(rRr))] );
hold on

plot(ts(1:length(rGr)),rGr,'g--')
plot(ts(1:length(rGr)),rGl,'g')
plot(ts(1:length(rGr)),rRr,'r--'); %vectorMag
[hAx,hLine1,hLine2] = plotyy(ts(1:length(rGl)),rRl,ts(1:length(rGl)),vR(1:length(rGl)));

ylabel(hAx(1),'Activity Vector Rate') % left y-axis 
ylabel(hAx(2),'Rotational Velocity') % right y-axis
xlim(hAx(1), [ ts(1) ts(end) ] )
xlim(hAx(2), [ ts(1) ts(end) ] )
ylim(hAx(1), [-ymax2 ymax2])
ylim(hAx(2), [-ymax ymax])
set(hLine1,'Color','r');
set(hLine2,'Color','k');
title('Rate of Change')
legend({strcat(green, ' R'), strcat(green, ' L'), strcat(red, 'R'), strcat(red, 'L'), 'vRot'});
set(hAx,{'ycolor'},{'k';'k'})

%%offset
subplot(2,1,2);
hold on

plot(ts, getCont(offsetr, 6), 'b--')
[hAx,hLine1,hLine2] = plotyy(ts,getCont(offsetl, 6),ts,vR);
ylim(hAx(1), [-4.5 4.5])
ylim(hAx(2), [-ymax ymax])
xlim(hAx(1), [ ts(1) ts(end) ] )
xlim(hAx(2), [ ts(1) ts(end) ] )
set(hAx,{'ycolor'},{'k';'k'})

ylabel(hAx(1),'Position Offset G-R') % left y-axis 
ylabel(hAx(2),'Rotational Velocity') % right y-axis
title('Offset')
set(hLine1,'Color','c');
set(hLine2,'Color','k');
legend({'Right', 'Left' 'vRot'});

%% plot pruned data

newact = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');


%rate of change
subplot(2,1,1);
ymax = max( max(vRp), abs(min(vRp)) )
ymax2 = max( [ max(abs(rGlp)) max(abs(rRlp)) max(abs(rGrp)) max(abs(rRrp))] )
hold on

plot(tp(1:length(rGrp)),rGrp,'g--');
plot(tp(1:length(rGrp)),rGlp,'g');
plot(tp(1:length(rGrp)),rRrp,'r--'); %rate of change

[hAx,hLine1,hLine2] = plotyy(tp(1:length(rRlp)),rRlp,tp(1:length(rRlp)),vRp(1:length(rRlp)));

ylabel(hAx(1),'Activity Vector Rate') % left y-axis 
ylabel(hAx(2),'Rotational Velocity') % right y-axis
xlim(hAx(1), [ tp(1) ts(end-1) ] )
xlim(hAx(2), [ tp(1) ts(end-1) ] )
ylim(hAx(1), [-ymax2 ymax2])
ylim(hAx(2), [-ymax ymax])
set(hLine1,'Color','r');
set(hLine2,'Color','k');
title('Rate of Change')
legend({strcat(green, ' R'), strcat(green, ' L'), strcat(red, 'R'), strcat(red, 'L'), 'vRot'});
set(hAx,{'ycolor'},{'k';'k'})

%%offset
subplot(2,1,2);
hold on
plot(tp(1:length(vRp)), getCont(offsetrp(1:length(vRp)), 6), '--b')
[hAx,hLine1,hLine2] = plotyy(tp(1:length(vRp)),getCont(offsetlp(1:length(vRp)), 6),tp(1:length(vRp)),vRp);
ylim(hAx(1), [-pi pi])
ylim(hAx(2), [-ymax ymax])
xlim(hAx(1), [ tp(1) tp(end) ] )
xlim(hAx(2), [ tp(1) tp(end) ] )
set(hAx,{'ycolor'},{'k';'k'})

ylabel(hAx(1),'Position Offset G-R') % left y-axis 
ylabel(hAx(2),'Rotational Velocity') % right y-axis
set(hLine1,'Color','c');
set(hLine2,'Color','k');
title('Offset')
legend({'Right', 'Left', 'vRot'});



%% Boxplots
%only do this for pruned data

xmin = min(vR);
xmax = max(vR);

newact2 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

subplot(2,2,1);

scatter(vRp, offsetrp(1:length(vRp)), 'b') %right offset
xlabel('vR')
ylabel('Offset G-R')

size(vRp)
size(offsetrp)

p1 = polyfit(vRp,offsetrp(1:length(vRp)),1);
yfit1 = polyval(p1,vRp);
yresid = offsetrp(1:length(vRp)) - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(offsetrp)-1) * var(offsetrp);
rsq1 = 1 - SSresid/SStotal;

slopes = [p1(1)];

hold on

scatter(vRp, offsetlp(1:length(vRp)), 'c') %left offset
xlabel('vR')

p2 = polyfit(vRp,offsetlp(1:length(vRp)),1);
yfit2 = polyval(p1,vRp);
yresid = offsetlp(1:length(vRp)) - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(offsetlp)-1) * var(offsetlp);
rsq2 = 1 - SSresid/SStotal;

slopes = [slopes p2(1)];

legend( {sprintf('Offset G-R R slope %.2f R2 %.2f', p1(1), rsq1),...
    sprintf('Offset G-R L slope %.2f R2 %.2f', p2(1), rsq2)} )

plot(vRp, yfit1, 'b')
plot(vRp, yfit2, 'c')
xlim([xmin xmax])


subplot(2,2,2) %rate of change of fluorescence direction right PB

scatter(vRp(1:length(rRrp)), rRrp, 'r')
hold on
scatter(vRp(1:length(rGrp)), rGrp, 'g')
xlabel('vR')
ylabel('Rate of change of fluorescence direction')

p1 = polyfit(vRp(1:length(rRrp)),rRrp,1);
yfit1 = polyval(p1,vRp(1:length(rRrp)));
yresid = rRrp - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(rRrp)-1) * var(rRrp);
rsq1 = 1 - SSresid/SStotal;

p2 = polyfit(vRp(1:length(rRrp)),rGrp,1);
yfit2 = polyval(p2,vRp(1:length(rRrp)));
yresid = rGrp - yfit2;
SSresid = sum(yresid.^2);
SStotal = (length(rGrp)-1) * var(rGrp);
rsq2 = 1 - SSresid/SStotal;

title('Rate of change Right')
legend( {sprintf('Red slope %.2f R2 %.2f', p1(1), rsq1), sprintf('Green slope %.2f R2 %.2f', p2(1), rsq2)} )

plot(vRp(1:length(rRrp)), yfit1, 'r')
plot(vRp(1:length(rRlp)), yfit2, 'g')
xlim([xmin xmax])

subplot(2,2,3) %rate of change of fluorescence direction left PB

scatter(vRp(1:length(rRlp)), rRlp, 'r')
hold on
scatter(vRp(1:length(rGlp)), rGlp, 'g')
xlabel('vR')
ylabel('Rate of change of fluorescence direction')

p1 = polyfit(vRp(1:length(rRlp)),rRlp,1);
yfit1 = polyval(p1,vRp(1:length(rRlp)));
yresid = rRlp - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(rRlp)-1) * var(rRlp);
rsq1 = 1 - SSresid/SStotal;

p2 = polyfit(vRp(1:length(rRlp)),rGlp,1);
yfit2 = polyval(p2,vRp(1:length(rRlp)));
yresid = rGlp - yfit2;
SSresid = sum(yresid.^2);
SStotal = (length(rGlp)-1) * var(rGlp);
rsq2 = 1 - SSresid/SStotal;

legend( {sprintf('Red slope %.2f R2 %.2f', p1(1), rsq1), sprintf('Green slope %.2f R2 %.2f', p2(1), rsq2)} )
title('Rate of change Left')
plot(vRp(1:length(rRrp)), yfit1, 'r')
plot(vRp(1:length(rRlp)), yfit2, 'g')
xlim([xmin xmax])

%% return data

cutoff = pi/3;

ominusR = mean(offsetrp( vRp < -cutoff )); %offset at +ve vs. -ve velocities
oplusR = mean(offsetrp( vRp > cutoff ));

ominusL = mean(offsetlp( vRp < -cutoff ));
oplusL = mean(offsetlp( vRp > cutoff ));

oAbsR = mean(abs(offsetrp)) %magnitude of offset
oAbsL = mean(abs(offsetlp))

actgR = mean( mGr( abs(vR) > cutoff )) / mean( mGr( abs(vR) < cutoff )) %intensity at high over low rotational velocities
actrR = mean( mRr( abs(vR) > cutoff )) / mean( mRr( abs(vR) < cutoff ))

actgL = mean( mGl( abs(vR) > cutoff )) / mean( mGl( abs(vR) < cutoff ))
actrL = mean( mRl( abs(vR) > cutoff )) / mean( mRl( abs(vR) < cutoff ))


%% save

name = strcat(dir, sprintf('fly%d _trial%d', fly, trial));

print(act, strcat(name, '_overview'), '-dpdf');
print(act2, strcat(name, '_apruned'), '-dpdf');
print(newact, strcat(name, '_pruned'),'-dpdf');
print(newact2, strcat(name, '_scatter'),'-dpdf');
