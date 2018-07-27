function [ominusR, oplusR, ominusL, oplusL, oAbsR, oAbsL, actgR, actrR, actgL, actrL, slopes] = plot_PB_Kris(dat, green, red, fly, trial, dir)

%%make a script for plot velocity, G and R as a function of t

%cond = FlyDatLoad(2)
%dat = cond{1}.allFlyData{1}.('Dark'){1}

vR = dat.positionDatMatch.vRot( dat.positionDatMatch.Closed(1:length(dat.positionDatMatch.vRot)) == 1 ) ; %CCW velocity
vF = dat.positionDatMatch.vF(dat.positionDatMatch.Closed(1:length(dat.positionDatMatch.vF))== 1); %forward velocity


ts = dat.positionDatMatch.OffsetRotMatch(:,1);
ts = ts(dat.positionDatMatch.Closed== 1);
%ts = ts(1:L); %timepoints for which we have v

L = length(dat.GROIaveMax);

heading = dat.positionDatMatch.PosRotMatch(1:L) * pi/180;
heading = heading(dat.positionDatMatch.Closed== 1)%CCW is positive

datG = flipud(dat.GROIaveMax)-1; %we're numbering from 1 at L9 to 18 at R9 on Tanya's diagram
datR = flipud(dat.RROIaveMax)-1; %this means CCW rotation (+ve vrot) gives higher numbers

size(datG);

s = size(datG);
s = s(1);

mG = []; %magnitude of vector sum of green EB activity
mR = [];

dG = []; %direction of vector sum
dR = [];

smooth = 3

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

L = length(datG);

[mGl, mGr, dGl, dGr] = getPBVec(datG)
[mRl, mRr, dRl, dRr] = getPBVec(datR)

mGt = mGr+mGl
mRt = mRr+mRl

L = length(vR); %number of data points we have velocity for

%%prune away stuff with too low intensity

prune = 0.20

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


rGr = getRates(dGr, ts, 9);%, tGs) %rate of change of direction of Green vector sum. +ve is moving to the right in the PB (e.g. CCW fly rotation)
rRr = getRates(dRr, ts, 9);%, tRs) %rate of change of direction of Red vector sum
rGrp = getRates(dGrp, tp, 9);
rRrp = getRates(dRrp, tp, 9);
rGl = getRates(dGl, ts, 9);
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

offsetr = getOffset(dGr, dRr, 9); %get offset of green relative to red
offsetrp = getOffset(dGrp, dRrp, 9);
offsetl = getOffset(dGl, dRl, 9);
offsetlp = getOffset(dGlp, dRlp, 9);

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

size(ts)
size(mGr)
size(mGl)
size(mRr)
size(mRl)

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

%we're looking at velocity so remove  bits where the fly is standing still

xmin = min(vR)
xmax = max(vR)

newact2 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

%include = [ abs(vR) > 0.1 ]
%vR = vR(include)
%offset = offset(include)
%rR = rR(include)
%rG = rG(include)
%{
subplot(2,2,1);

scatter(vR, offset)
xlabel('vRot')
ylabel('Offset G-R')

p1 = polyfit(vR,offset,1)
yfit1 = polyval(p1,vR);
yresid = offset - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(offset)-1) * var(offset);
rsq1 = 1 - SSresid/SStotal

legend( {sprintf('Offset G-R slope %.2f R2 %.2f', p1(1), rsq1)} )
hold on
plot(vR, yfit1)
xlim([xmin xmax])
%}

subplot(2,2,1);

%include = (abs(vRp) > 0.1)
%vRp = vRp(include)
%offsetp = offsetp(include)
%rRp = rRp(include)
%rGp = rGp(include)

scatter(vRp, offsetrp(1:length(vRp)), 'b')
xlabel('vR')
ylabel('Offset G-R')

p1 = polyfit(vRp,offsetrp(1:length(vRp)),1)
yfit1 = polyval(p1,vRp);
yresid = offsetrp(1:length(vRp)) - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(offsetrp)-1) * var(offsetrp);
rsq1 = 1 - SSresid/SStotal

slopes = [p1(1)]

hold on


scatter(vRp, offsetlp(1:length(vRp)), 'c')
xlabel('vR')

p2 = polyfit(vRp,offsetlp(1:length(vRp)),1)
yfit2 = polyval(p1,vRp);
yresid = offsetlp(1:length(vRp)) - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(offsetlp)-1) * var(offsetlp);
rsq2 = 1 - SSresid/SStotal

slopes = [slopes p2(1)]

legend( {sprintf('Offset G-R R slope %.2f R2 %.2f', p1(1), rsq1),...
    sprintf('Offset G-R L slope %.2f R2 %.2f', p2(1), rsq2)} )

plot(vRp, yfit1, 'b')
plot(vRp, yfit2, 'c')
xlim([xmin xmax])

%{
subplot(2,2,3)

scatter(vR(1:length(rR)), rRr, 'r')
hold on
scatter(vR(1:length(rG)), rGr, 'g')

legend( {'EPG', 'PEN1'} )
xlabel('vR')
ylabel('Rate of change of fluorescence direction')

p1 = polyfit(vR(1:length(rR)),rR,1);
yfit1 = polyval(p1,vR(1:length(rR)));
yresid = rR - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(rR)-1) * var(rR);
rsq1 = 1 - SSresid/SStotal;

scatter(vR(1:length(rR)), rR, 'r')
hold on
scatter(vR(1:length(rG)), rG, 'g')

legend( {'EPG', 'PEN1'} )
xlabel('vR')
ylabel('Rate of change of fluorescence direction')

p2 = polyfit(vR(1:length(rR)),rG,1);
yfit2 = polyval(p2,vR(1:length(rR)));
yresid = rG - yfit2;
SSresid = sum(yresid.^2);
SStotal = (length(rG)-1) * var(rG);
rsq2 = 1 - SSresid/SStotal;

legend( {sprintf('Red slope %.2f R2 %.2f', p1(1), rsq1), sprintf('Green slope %.2f R2 %.2f', p2(1), rsq2)} )
hold on
plot(vR(1:length(rR)), yfit1, 'r')
plot(vR(1:length(rR)), yfit2, 'g')
xlim([xmin xmax])
%}

subplot(2,2,2)


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



subplot(2,2,3)


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

cutoff = pi/3

ominusR = mean(offsetrp( vRp < -cutoff ))
oplusR = mean(offsetrp( vRp > cutoff ))

ominusL = mean(offsetlp( vRp < -cutoff ))
oplusL = mean(offsetlp( vRp > cutoff ))

oAbsR = mean(abs(offsetrp))
oAbsL = mean(abs(offsetlp))

actgR = mean( mGr( abs(vR) > cutoff )) / mean( mGr( abs(vR) < cutoff ))
actrR = mean( mRr( abs(vR) > cutoff )) / mean( mRr( abs(vR) < cutoff ))

actgL = mean( mGl( abs(vR) > cutoff )) / mean( mGl( abs(vR) < cutoff ))
actrL = mean( mRl( abs(vR) > cutoff )) / mean( mRl( abs(vR) < cutoff ))


%% save

name = strcat(dir, sprintf('fly%d _trial%d', fly, trial));

print(act, strcat(name, '_overview'), '-dpdf');
print(act2, strcat(name, '_apruned'), '-dpdf');
print(newact, strcat(name, '_pruned'),'-dpdf');
print(newact2, strcat(name, '_scatter'),'-dpdf');
