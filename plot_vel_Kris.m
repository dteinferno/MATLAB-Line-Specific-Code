function [mminus, stdminus, mplus, stdplus, actg, actr, slope] = plot_vel_Kris(dat, green, red, fly, trial, dir)

%%make a script for plot velocity, G and R as a function of t

%cond = FlyDatLoad(2)
%dat = cond{1}.allFlyData{1}.('Dark'){1}

vR = dat.positionDatMatch.vRot( dat.positionDatMatch.Closed(1:length(dat.positionDatMatch.vRot)) == 1 ) ; %rotational velocity
vF = dat.positionDatMatch.vF(dat.positionDatMatch.Closed(1:length(dat.positionDatMatch.vF))== 1); %forward velocity


ts = dat.positionDatMatch.OffsetRotMatch(:,1);
ts = ts(dat.positionDatMatch.Closed== 1);
%ts = ts(1:L); %timepoints for which we have v

L = length(dat.GROIaveMax);

heading = dat.positionDatMatch.PosRotMatch(1:L) * pi/180;
heading = heading(dat.positionDatMatch.Closed== 1)%fly heading actually not heading but expected bump position?

datG = dat.GROIaveMax;
datR = dat.RROIaveMax;

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
    newR = [];
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

for i = 1:L;
    
    [dGi, mGi] = getVecSum( datG(:, i) );
    dG = vertcat(dG, dGi);
    mG = vertcat(mG, mGi);
    
    [dRi, mRi] = getVecSum( datR(:, i) );
    dR = vertcat(dR, dRi);
    mR = vertcat(mR, mRi);
    
end
dR = dR-pi;
dG = dG-pi;

L = length(vR); %number of data points we have velocity for

%%prune away stuff with too low intensity

prune = 0.35

if prune > 0
    rlim = prune*mean(mR);
    glim = prune*mean(mG);
    
    dGp = dG( mR > rlim );
    tp = ts( mR > rlim );
    vRp = vR( mR(1:end-1) > rlim ); 
    mGp = mG( mR > rlim );
    dRp = dR( mR > rlim );
    mRp = mR( mR > rlim );
    dGp = dGp( mGp > glim );
    tp = tp( mGp > glim );
    vRp = vRp( mGp(1:end-1) > glim ); 
    dRp = dRp( mGp > glim );
    mRp = mRp( mGp > glim );
    mGp = mGp( mGp > glim );
end


rG = getRates(dG, ts, 2*pi);%, tGs) %rate of change of direction of Green vector sum
rR = getRates(dR, ts, 2*pi);%, tRs) %rate of change of direction of Red vector sum
rGp = getRates(dGp, tp, 2*pi);
rRp = getRates(dRp, tp, 2*pi);

%cut off last value to amke copatible with our rate data
dG = dG(1:L);
dR = dR(1:L);
mG = mG(1:L);
mR = mR(1:L);
ts = ts(1:L);
heading = heading(1:L);

offset = getOffset(dG, dR, 2*pi); %get offset of green relative to red
offsetp = getOffset(dGp, dRp, 2*pi);

%% plot raw data

% Plot the profiles
act = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

% Plot the activity
subplot(3,1,1);
hold on;
plot(ts,getCont(heading, pi),'b'); %time vs heading
plot(ts,getCont(dG, 1.2*pi),'g'); %vectorPos
plot(ts,getCont(dR, 1.2*pi),'r');%vectorPoos
legend({'heading',green, red});
ylabel('position (rad)');
ylim([-pi pi] );
xlim( [ ts(1) ts(end) ] );

%%forward velocity
subplot(3,1,2);
absvF = abs(vF);
hold on
plot(ts,mG,'g'); %vectorMag

[hAx,hLine1,hLine2] = plotyy(ts,mR,ts,absvF);
set(hAx,{'ycolor'},{'k';'k'})

ylabel(hAx(1),'Position Vector Magnitude') % left y-axis 
ylabel(hAx(2),'Forward Velocity Magnitude') % right y-axis
xlim(hAx(1), [ ts(1) ts(end) ] )
xlim(hAx(2), [ ts(1) ts(end) ] )
ylim(hAx(1), [0 max( max(mR), max(mG) )])

min(absvF)
max(absvF)
vF

ylim(hAx(2), [min( absvF ) max(absvF)])
set(hLine1,'Color','r');
set(hLine2,'Color','k');
legend({green, red, 'vF'});

%%rotatinal velocity
subplot(3,1,3);
absvR = abs(vR);
hold on
plot(ts,mG,'g'); %vectorMag
[hAx,hLine1,hLine2] = plotyy(ts,mR,ts,absvR);
set(hAx,{'ycolor'},{'k';'k'})

ylabel(hAx(1),'Activity Vector Magnitude') % left y-axis 
ylabel(hAx(2),'Rotational Velocity Magnitude') % right y-axis
xlim(hAx(1), [ ts(1) ts(end) ] )
xlim(hAx(2), [ ts(1) ts(end) ] )
ylim(hAx(1), [0 max( max(mR), max(mG) )])
ylim(hAx(2), [min( absvR ) max(absvR)])
set(hLine1,'Color','r');
set(hLine2,'Color','k');
legend({green, red, 'vR'});

%% plot vs offset and velocity
act2 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

%rate of change
subplot(2,1,1);
ymax = max( max(vR), abs(min(vR)) );
ymax2 = max( max(abs(rG)), max(abs(rR)));
hold on

plot(ts(1:length(rG)),rG,'g'); %vectorMag
[hAx,hLine1,hLine2] = plotyy(ts(1:length(rG)),rR,ts(1:length(rG)),vR(1:length(rG)));

ylabel(hAx(1),'Activity Vector Rate of Change') % left y-axis 
ylabel(hAx(2),'Rotational Velocity') % right y-axis
xlim(hAx(1), [ ts(1) ts(end) ] )
xlim(hAx(2), [ ts(1) ts(end) ] )
ylim(hAx(1), [-ymax2 ymax2])
ylim(hAx(2), [-ymax ymax])
set(hLine1,'Color','r');
set(hLine2,'Color','k');
legend({green, red, 'vR'});
set(hAx,{'ycolor'},{'k';'k'})

%%offset
subplot(2,1,2);

[hAx,hLine1,hLine2] = plotyy(ts,offset,ts,vR);
ylim(hAx(1), [-pi pi])
ylim(hAx(2), [-ymax ymax])
xlim(hAx(1), [ ts(1) ts(end) ] )
xlim(hAx(2), [ ts(1) ts(end) ] )
set(hAx,{'ycolor'},{'k';'k'})

ylabel(hAx(1),'Position Vector Offset G-R') % left y-axis 
ylabel(hAx(2),'Rotational Velocity') % right y-axis
set(hLine1,'Color','r');
set(hLine2,'Color','k');
legend({'Offset G-R', 'vR'});

%% plot pruned data

newact = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');


%rate of change
subplot(2,1,1);
ymax = max( max(vRp), abs(min(vRp)) )
ymax2 = max( max(abs(rGp)), max(abs(rRp)))
hold on
plot(tp(1:length(rGp)),rGp,'g'); %rate of change

[hAx,hLine1,hLine2] = plotyy(tp(1:length(rRp)),rRp,tp(1:length(rRp)),vRp(1:length(rRp)));

ylabel(hAx(1),'Activity Vector Rate of Change') % left y-axis 
ylabel(hAx(2),'Rotational Velocity') % right y-axis
xlim(hAx(1), [ tp(1) ts(end-1) ] )
xlim(hAx(2), [ tp(1) ts(end-1) ] )
ylim(hAx(1), [-ymax2 ymax2])
ylim(hAx(2), [-ymax ymax])
set(hLine1,'Color','r');
set(hLine2,'Color','k');
legend({green, red, 'vR'});
set(hAx,{'ycolor'},{'k';'k'})

%%offset
subplot(2,1,2);

[hAx,hLine1,hLine2] = plotyy(tp(1:length(vRp)),offsetp(1:length(vRp)),tp(1:length(vRp)),vRp);
ylim(hAx(1), [-pi pi])
ylim(hAx(2), [-ymax ymax])
xlim(hAx(1), [ tp(1) tp(end) ] )
xlim(hAx(2), [ tp(1) tp(end) ] )
set(hAx,{'ycolor'},{'k';'k'})

ylabel(hAx(1),'Position Vector Offset G-R') % left y-axis 
ylabel(hAx(2),'Rotational Velocity') % right y-axis
set(hLine1,'Color','r');
set(hLine2,'Color','k');
legend({'Offset G-R', 'vR'});

%% Boxplots

%we're looking at velocity so remove  bits where the fly is standing still

xmin = min(vR);
xmax = max(vR);

newact2 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

%include = [ abs(vR) > 0.1 ]
%vR = vR(include)
%offset = offset(include)
%rR = rR(include)
%rG = rG(include)

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

subplot(2,2,2);

%include = (abs(vRp) > 0.1)
%vRp = vRp(include)
%offsetp = offsetp(include)
%rRp = rRp(include)
%rGp = rGp(include)

scatter(vRp, offsetp(1:length(vRp)))
xlabel('vR')
ylabel('Offset G-R')

p1 = polyfit(vRp,offsetp(1:length(vRp)),1)
yfit1 = polyval(p1,vRp);
yresid = offsetp(1:length(vRp)) - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(offsetp)-1) * var(offsetp);
rsq1 = 1 - SSresid/SStotal

slope = p1(1)

legend( {sprintf('Offset G-R slope %.2f R2 %.2f', p1(1), rsq1)} )
hold on
plot(vRp, yfit1)
xlim([xmin xmax])

subplot(2,2,3)

scatter(vR(1:length(rR)), rR, 'r')
hold on
scatter(vR(1:length(rG)), rG, 'g')

legend( {'EPG', 'PEN1'} )
xlabel('vR')
ylabel('Rate of change of fluorescence direction')

p1 = polyfit(vR(1:length(rR)),rR,1);
yfit1 = polyval(p1,vR(1:length(rR)));
yresid = rR - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(rR)-1) * var(rR);
rsq1 = 1 - SSresid/SStotal;

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

subplot(2,2,4)


scatter(vRp(1:length(rRp)), rRp, 'r')
hold on
scatter(vRp(1:length(rGp)), rGp, 'g')
xlabel('vR')
ylabel('Rate of change of fluorescence direction')

p1 = polyfit(vRp(1:length(rRp)),rRp,1);
yfit1 = polyval(p1,vRp(1:length(rRp)));
yresid = rRp - yfit1;
SSresid = sum(yresid.^2);
SStotal = (length(rRp)-1) * var(rRp);
rsq1 = 1 - SSresid/SStotal;

p2 = polyfit(vRp(1:length(rRp)),rGp,1);
yfit2 = polyval(p2,vRp(1:length(rRp)));
yresid = rGp - yfit2;
SSresid = sum(yresid.^2);
SStotal = (length(rGp)-1) * var(rGp);
rsq2 = 1 - SSresid/SStotal;

legend( {sprintf('Red slope %.2f R2 %.2f', p1(1), rsq1), sprintf('Green slope %.2f R2 %.2f', p2(1), rsq2)} )
hold on
plot(vRp(1:length(rRp)), yfit1, 'r')
plot(vRp(1:length(rRp)), yfit2, 'g')
xlim([xmin xmax])



%% return data

cutoff = pi/3

mminus = mean(offsetp( vRp < -cutoff ))
stdminus = std(offsetp( vRp < -cutoff ))
mplus = mean(offsetp( vRp > cutoff ))
stdplus = std(offsetp( vRp > cutoff ))

actg = mean( mG( abs(vR) > cutoff )) / mean( mG( abs(vR) < cutoff ))
actr = mean( mR( abs(vR) > cutoff )) / mean( mR( abs(vR) < cutoff ))


%% save

name = strcat(dir, sprintf('fly%d _trial%d', fly, trial));

print(act, strcat(name, '_overview'), '-dpdf');
print(act2, strcat(name, '_apruned'), '-dpdf');
print(newact, strcat(name, '_pruned'),'-dpdf');
print(newact2, strcat(name, '_scatter'),'-dpdf');



