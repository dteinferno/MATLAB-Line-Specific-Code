function [taus, tPts] = getDecayDelayPB(dir, cond, l, delay, green, nWedge)
%Plots the decay of fluorescent activity when the fly stops moving and the
%re-emergence of activity when the fly starts moving again. Fits an
%exponential to the decay and sigmoidal to the re-emergence.
%type is either max(use max glomerulus) or mean


%%

%l is length of stretch that must be zero
%delay is length of stretch to consider for re-emergence of fluorescence

%find periods of 5+ time points of low vRot+vF
%align beginning and look at fluorescence. Normalize to beginning
%align end and look at fluorescence. normalize to end

threshR = 0; %want the fly to actually stand still
threshF = 0; %1/100;
pre1 = 2; %dictates how many prior points we include.
pre2 = 1;


name = cond{1}.name;

pers = { 'dark', 'OL', 'CL' };

taus = containers.Map(pers, {[], [], []});

for per = 1:3

    [DFs, data, flyinds] = getData(cond, pers{per}, 7, 0, nWedge); %smooth fluorescence but not vel. Makes it somewhat less noisy

    DFsL = DFs(1:9, :);
    DFsR = DFs(10:18, :);
    
    maxL = max(DFsL, [], 1);
    maxR = max(DFsR, [], 1);
    
    intL = mean( maxk(DFsL, nWedge, 1), 1);
    intR = mean( maxk(DFsR, nWedge, 1), 1);
    
    tPts = data{1};
    vR = abs( data{2} );
    vF = abs( data{3} );
    int = data{4};
    mags = data{5};
    dirs = data{6};
    vDF = abs(getRates(dirs, tPts, 2*pi) ); %get rates

    decaysL = [];
    delaysL = [];
    decaysR = [];
    delaysR = [];
    vdecs = [];
    vdels = [];
    vFdecs = [];
    vFdels = [];%magnitude of rate of change of direction of PVA
    magdecs = [];%magnitude of PVA
    magdels = [];
    
    if length(tPts) > 6
        xvals1 = [0:l+pre1-1].*(tPts(6)-tPts(5)); %length of decays
        xvals2 = [0:delay+pre2-1].*(tPts(6)-tPts(5)); %length of delays
    end

    prevind = 1;
    
    for i = 3:length(vR)-(l-1)
        %need continuous stationary stretch of length l starting at i
        if any(flyinds == i) || i == length(vR)-(l-1) %new trial; add taus from previous trial
            display('starting new fly')
            s = size(decaysL);
            if s(2) >= prevind %have added bouts of standing
                mdecL = mean(decaysL(:, prevind:s(2)), 2); %mean intensity
                mdecR = mean(decaysR(:, prevind:s(2)), 2);
                fdecL = fit(xvals1', mdecL, 'exp1');
                fdecR = fit(xvals1', mdecR, 'exp1');
                
                taus(pers{per}) = [taus(pers{per}); -fdecL.b -fdecR.b]; %rate constants

                prevind = s(2)+1;
            end
        end
        if all(vR(i:i+(l-1)) <= threshR) && all(vF(i:i+(l-1)) <= threshF ) && ...
                all(vR(i-pre1:i-1) > threshR) && all(vF(i-pre1:i-1) > threshF ) && ...
                tPts(i+(l-1)) > tPts(i-pre1)
            
            decaysL = [decaysL transpose( intL(i-pre1:i+(l-1))/intL(i-pre1) ) ]; %normalize to first data point
            decaysR = [decaysR transpose( intR(i-pre1:i+(l-1))/intR(i-pre1) ) ];
            
            %{
            if strcmp(method, 'mean')
                decaysL = [decaysL transpose( intL(i-pre1:i+(l-1))/intL(i-pre1) ) ]; %normalize to first data point
                decaysR = [decaysR transpose( intR(i-pre1:i+(l-1))/intR(i-pre1) ) ];
            elseif strcmp(method, 'max')
                decaysL = [decaysL transpose( maxL(i-pre1:i+(l-1))/maxL(i-pre1) ) ]; %normalize to first data point
                decaysR = [decaysR transpose( maxR(i-pre1:i+(l-1))/maxR(i-pre1) ) ];
            else
                display('method should be mean or max')
            end
            %}
                
            vFdecs = [vFdecs vDF(i-pre1:i+(l-1)) ]; %no normaliztion
            magdecs = [magdecs transpose( mags(i-pre1:i+(l-1)) ) ];%/ mags(i-pre1) ]; %normalize
            vdecs = [vdecs transpose( vR(i-pre1:i+(l-1)) ) ];

            j = l;
            while vR(i+j) <= threshR && vF(i+j) <= threshF %find the end of the stationary stretch
                j = j+1;
            end

            if length(vR) > i+j+(delay-2) %don't want to reach end of list
                %display('adding delay')
                %vR( i+j : i+j+(delay-1) )
                %vF( i+j : i+j+(delay-1) )

                %want continuous activity
                if all( vR( i+j : i+j+(delay-1) ) > threshR ) && all( vF( i+j : i+j+(delay-1) ) > threshF )  && ...
                        tPts(i+j+(delay-1)) > tPts(i+j-pre2)
                    %display('still adding delay')
                    
                    delaysL = [delaysL transpose( intL( i+j-pre2 : i+j+(delay-1) ) / intL(i+j-pre2) )];
                    delaysR = [delaysR transpose( intR( i+j-pre2 : i+j+(delay-1) ) / intR(i+j-pre2) )];
                    
                    %{
                    if strcmp(method, 'mean')
                        delaysL = [delaysL transpose( intL( i+j-pre2 : i+j+(delay-1) ) / intL(i+j-pre2) )];
                        delaysR = [delaysR transpose( intR( i+j-pre2 : i+j+(delay-1) ) / intR(i+j-pre2) )];
                    elseif strcmp(method, 'max')
                        delaysL = [delaysL transpose( maxL( i+j-pre2 : i+j+(delay-1) ) / maxL(i+j-pre2) )];
                        delaysR = [delaysR transpose( maxR( i+j-pre2 : i+j+(delay-1) ) / maxR(i+j-pre2) )];
                    else
                        display('method should be mean or max')
                    end
                    %}
                    
                    vFdels = [vFdels vDF( i+j-pre2 : i+j+(delay-1) ) ];
                    magdels = [magdels transpose( mags( i+j-pre2 : i+j+(delay-1) ))];% / mags(i+j-pre2) )];
                    vdels = [vdels transpose( vR( i+j-pre2 : i+j+(delay-1) ) ) ];
                end

            end                
        end
    end

    if length(vdecs) > 0
    
        %vdecs
        %decaysL
        %ecaysR
        %vFdecs
        %magdecs

        %vdels
        %delaysL
        %delaysR
        %Fdels
        %magdels

        %fit delay to sigmoidal function of form
        %f(x) = 1 + a / ( 1 + exp[-b(x-c)] )
        %delay is given by c
        %max slope if given by b*c

        %fot
        fo = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0,0],...
            'Upper',[5,100,1],...
            'StartPoint',[0.5 4 0.5]);
        ft = fittype('y0+a/(1+exp(-b*(x-c)))','options',fo,'problem','y0'); %allow us to specify y0 as first value of list

        mvdec = mean(vdecs, 2); %mean velocity
        mdecL = mean(decaysL,2); %mean intensity
        mdecR = mean(decaysR,2); %mean intensity
        
        vFdecs(vFdecs == Inf) = mean(vFdecs(isfinite(vFdecs))); %deal with infinite values...
        mvFdec = mean(vFdecs, 2); %mean rate of change of fluorescence direction
        mmagdec = mean(magdecs, 2); %mean of |PVA|

        mvdel = mean(vdels, 2);
        mdelL = mean(delaysL,2);
        mdelR = mean(delaysR,2);
        
        mvFdel = mean(vFdels, 2);
        mmagdel = mean(magdels, 2);

        newfig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

        fdecL = fit(xvals1', mdecL, 'exp1');
        [fdelL, gofdelL] = fit(xvals2', mdelL, ft, 'problem', mdelL(1));     
        fitdecL = fdecL.a * exp( fdecL.b * xvals1); %we fit to exp(kx) NOT exp(-kx)
        fitdelL = mdelL(1) + fdelL.a * ( 1 + exp( -fdelL.b * (xvals2 - fdelL.c) ) ).^(-1); %fit sigmoidal
        
        fdecR = fit(xvals1', mdecR, 'exp1');
        [fdelR, gofdelR] = fit(xvals2', mdelR, ft, 'problem', mdelR(1));     
        fitdecR = fdecR.a * exp( fdecR.b * xvals1); %we fit to exp(kx) NOT exp(-kx)
        fitdelR = mdelR(1) + fdelR.a * ( 1 + exp( -fdelR.b * (xvals2 - fdelR.c) ) ).^(-1); %fit sigmoidal


        %plot intensities decay
        subplot(3,2,1)
        hold on

        p1 = plot(xvals1, mdecL, 'k', 'linewidth', 5);
        title( strcat('intensity-decays-',green, '-', name,'-',pers(per)) )
        p2 = plot(xvals1, fitdecL, 'b--', 'linewidth', 3);
        legend( { 'Mean', strcat('tau =',num2str(-1/fdecL.b, 3)) } )
        plot(xvals1, decaysL)
        uistack([p2 p1], 'top')
        xlim( [0 max(xvals1)] )
        ylim([0 2])

        %plot intensities delay
        subplot(3,2,2)
        hold on
        title( strcat('intensity-delays-',green, '-', name,'-',pers(per)) )
        p1 = plot(xvals2, mdelL, 'k', 'linewidth', 5);
        p2 = plot(xvals2, fitdelL, 'b--', 'linewidth', 3);
        legend( { 'Mean', strcat('delay =',num2str(fdelL.c, 3),' slope=',num2str(0.25*fdelL.a * fdelL.b, 3))    })%,...
            %[num2str(fdel.a) '-' num2str(fdel.b) '-' num2str(fdel.c)] } )
        plot(xvals2, delaysL)
        uistack([p2 p1], 'top')
        xlim( [0 max(xvals2)] )
        ylim([0 5])
        
        %plot intensities decay
        subplot(3,2,3)
        hold on

        p1 = plot(xvals1, mdecR, 'k', 'linewidth', 5);
        title( strcat('intensity-decays-',green, '-', name,'-',pers(per)) )
        p2 = plot(xvals1, fitdecR, 'b--', 'linewidth', 3);
        legend( { 'Mean', strcat('tau =',num2str(-1/fdecR.b, 3)) } )
        plot(xvals1, decaysR)
        uistack([p2 p1], 'top')
        xlim( [0 max(xvals1)] )
        ylim([0 2])

        %plot intensities delay
        subplot(3,2,4)
        hold on
        title( strcat('intensity-delays-',green, '-', name,'-',pers(per)) )
        p1 = plot(xvals2, mdelR, 'k', 'linewidth', 5);
        p2 = plot(xvals2, fitdelR, 'b--', 'linewidth', 3);
        legend( { 'Mean', strcat('delay =',num2str(fdelR.c, 3),' slope=',num2str(0.25*fdelR.a * fdelR.b, 3))    })%,...
            %[num2str(fdel.a) '-' num2str(fdel.b) '-' num2str(fdel.c)] } )
        plot(xvals2, delaysR)
        uistack([p2 p1], 'top')
        xlim( [0 max(xvals2)] )
        ylim([0 5])





        newfig.PaperUnits = 'inches';
        newfig.PaperPosition = [0 0 8 11.5];

        %dest = strcat( dir, 'decay_delay_', name, '_', pers(per), '_', conds{condid} )
        dest = strcat(dir, 'decayDelay_',green, '_', name,'_',pers{per},'_','nWedge',num2str(nWedge));

        print(newfig, dest, '-dpdf');
    end


end