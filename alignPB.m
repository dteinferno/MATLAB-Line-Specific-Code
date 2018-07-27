function [mGL, stdGL, mRL, stdRL, mGR, stdGR, mRR, stdRR, dirGL, magGL, dirRL, magRL, dirGR, magGR,...
            dirRR, magRR] = alignPB(dat1, dat2, n1, n2)

%datG = flipud(dat.GROIaveMax); %we're numbering from 1 at L9 to 18 at R9 on Tanya's diagram
%datR = flipud(dat.RROIaveMax); %this means CCW rotation (+ve vrot) gives higher numbers

GalR = [];
GalL = [];
RalR = [];
RalL = [];

PEN = {'PEN1', 'PEN2'}
EPG = {'EPG', 'PEG'}


s = size(dat1);

for i = 1:s(2);
    
    l1 = 8
    l2 = 8
    
    if strcmp(n1, 'PEN1') || strcmp(n1, 'PEN2')
        GL = dat1(1:8,i); %outer glomeruli
        GR = dat1(11:18,i);
        sL=5;
        sR=4;
    elseif strcmp(n1, 'EPG') || strcmp(n1, 'PEG')
        GL = dat1(2:9,i); %inner glomeruli
        GR = dat1(10:17,i);
        sL=4;
        sR=5;
    else
        GL = dat1(1:9,i); %all glomeruli
        GR = dat1(10:18,i);
        l1 = 9; %need to rotate by 9 if D7
        sL=5;
        sR=5;
    end

    if strcmp(n2, 'PEN1') || strcmp(n2, 'PEN2')
        RL = dat2(1:8,i); %outer glomeruli
        RR = dat2(11:18,i);
    elseif strcmp(n2, 'EPG') || strcmp(n2, 'PEG')
        RL = dat2(2:9,i); %inner glomeruli
        RR = dat2(10:17,i);
    else
        RL = dat2(1:9,i); %all glomeruli
        RR = dat2(10:18,i);
        l2 = 9; %need to rotate by 9 if D7
    end
        
    %% Shift stuff
    newGL = GL;
    newRL = RL;
    newGR = GR;
    newRR = RR;
    
    [valL idxL] = max(GL); %index in truncated list
    shiftL = sL-idxL;
    
    [valR idxR] = max(GR); %index in truncated list
    shiftR = sR-idxR;
    
    for j = 1:8
        kL = j+shiftL;
        kR = j+shiftR;
        
        if kL > l1
            newGL(kL-l1) = GL(j);     
        elseif kL < 1
            newGL(kL+l1) = GL(j);
        else    
            newGL(kL) = GL(j);
        end
        if kL > l2
            newRL(kL-l2) = RL(j);      
        elseif kL < 1
            newRL(kL+l2) = RL(j);
        else    
            newRL(kL) = RL(j);
        end
        
        if kR > l1
            newGR(kR-l1) = GR(j);
        elseif kR < 1
            newGR(kR+l1) = GR(j);
        else
            newGR(kR) = GR(j);
        end
        if kR > l2
            newRR(kR-l2) = RR(j);
        elseif kR < 1
            newRR(kR+l2) = RR(j);
        else
            newRR(kR) = RR(j);
        end
    end

    %% need to fix element 9 if D7
    if strcmp(n1, 'D7')
        kL = 9+shiftL;
        kR = 9+shiftR;
        
        if kL > l1
            newGL(kL-l1) = GL(9);     
        elseif kL < 1
            newGL(kL+l1) = GL(9);
        else    
            newGL(kL) = GL(9);
        end
        if kR > l1
            newGR(kR-l1) = GR(j);
        elseif kR < 1
            newGR(kR+l1) = GR(j);
        else
            newGR(kR) = GR(j);
        end   
    end
    if strcmp(n2, 'D7')
        kL = 9+shiftL;
        kR = 9+shiftR;
        if kL > l2
            newRL(kL-l2) = RL(j);      
        elseif kL < 1
            newRL(kL+l2) = RL(j);
        else    
            newRL(kL) = RL(j);
        end
        if kR > l2
            newRR(kR-l2) = RR(j);
        elseif kR < 1
            newRR(kR+l2) = RR(j);
        else
            newRR(kR) = RR(j);
        end
        
    end
    %% Add data
    if strcmp(n1, 'PEN1') || strcmp(n1, 'PEN2')
        GalL(:,i) = vertcat(newGL, 0);
        GalR(:,i) = vertcat(0, newGR);
    elseif strcmp(n1, 'EPG') || strcmp(n1, 'PEG')
        GalL(:,i) = vertcat(0, newGL);
        GalR(:,i) = vertcat(newGR, 0);
    else
        GalL(:,i) = newGL;
        GalR(:,i) = newGR;
    end
        
    if strcmp(n2, 'PEN1') || strcmp(n2, 'PEN2')
        RalL(:,i) = vertcat(newRL, 0);
        RalR(:,i) = vertcat(0, newRR);
    elseif strcmp(n2, 'EPG') || strcmp(n2, 'PEG')
        RalL(:,i) = vertcat(0, newRL);
        RalR(:,i) = vertcat(newRR, 0);
    else
        RalL(:,i) = newRL;
        RalR(:,i) = newRR;
    end
end

%% Get mean values
mGL = [];
stdGL = [];
mRL = [];
stdRL = [];
mGR = [];
stdGR = [];
mRR = [];
stdRR = [];

for i = 1:9
    mGL = [mGL mean(GalL(i,:))];
    stdGL = [stdGL std(GalL(i,:))];
    mRL = [mRL mean(RalL(i,:))];
    stdRL = [stdRL std(RalL(i,:))];
    mGR = [mGR mean(GalR(i,:))];
    stdGR = [stdGR std(GalR(i,:))];
    mRR = [mRR mean(RalR(i,:))];
    stdRR = [stdRR std(RalR(i,:))];
end

if strcmp(n1, 'PEN1') || strcmp(n1, 'PEN2')
    [dirGL, magGL] = getVecSum(mGL(1:8));
    [dirGR, magGR] = getVecSum(mGR(2:9));
    dirGL = dirGL * 8/(2*pi);
    dirGR = dirGR * 8/(2*pi)+1;%need to add in the discounted glomerulus

elseif strcmp(n1, 'EPG') || strcmp(n1, 'PEG')
    [dirGL, magGL] = getVecSum(mGL(2:9));
    [dirGR, magGR] = getVecSum(mGR(1:8));
    dirGL = dirGL * 8/(2*pi)+1;%need to add in the discounted glomerulus
    dirGR = dirGR * 8/(2*pi);
else
    [dirGL, magGL] = getVecSum(mGL);
    [dirGR, magGR] = getVecSum(mGR);
    dirGL = dirGL * 9/(2*pi);%D7 we consider intact
    dirGR = dirGR * 9/(2*pi);
end

if strcmp(n2, 'PEN1') || strcmp(n2, 'PEN2') %get vector sum from circle of 8 pieces, then add missing glomerulus
    %effectively we imagine that the PB is 2x8 glomeruli and then just add
    %the offset at the end
    [dirRL, magRL] = getVecSum(mRL(1:8));
    [dirRR, magRR] = getVecSum(mRR(2:9));
    dirRL = dirRL * 8/(2*pi);
    dirRR = dirRR * 8/(2*pi)+1;%need to add in the discounted glomerulus

elseif strcmp(n2, 'EPG') || strcmp(n2, 'PEG')
    [dirRL, magRL] = getVecSum(mRL(2:9));
    [dirRR, magRR] = getVecSum(mRR(1:8));
    dirRL = dirRL * 8/(2*pi)+1;%need to add in the discounted glomerulus
    dirRR = dirRR * 8/(2*pi);
else
    [dirRL, magRL] = getVecSum(mRL);
    [dirRR, magRR] = getVecSum(mRR);
    dirRL = dirRL * 9/(2*pi);%D7 we consider intact
    dirRR = dirRR * 9/(2*pi);
end

magRL = magRL / sum(mRL) * 3; %get 25 times PVA
magGL = magGL / sum(mGL) * 3;
magRR = magRR / sum(mRR) * 3; %get 25 times PVA
magGR = magGR / sum(mGR) * 3;


%{
fig = figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,1,1)
hold on
plot(1:9, mGL, 'g')
plot(1:9, mRL, 'r')

plot(1:9, mGL+stdGL, 'g--')
plot(1:9, mGL-stdGL, 'g--')

plot(1:9, mRL+stdRL, '--r')
plot(1:9, mRL-stdRL, '--r')

xlim([1 9])
title('Left PB')
legend({n1, n2})

subplot(2,1,2)
hold on

plot(1:9, mGR, 'g')
plot(1:9, mRR, 'r')

plot(1:9, mGR+stdGR, 'g--')
plot(1:9, mGR-stdGR, 'g--')

plot(1:9, mRR+stdRR, '--r')
plot(1:9, mRR-stdRR, '--r')

xlim([1 9])
title('Right PB')
legend({n1, n2})

%}
