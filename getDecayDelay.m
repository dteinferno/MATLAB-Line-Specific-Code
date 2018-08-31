function taus = getDecayDelay(dir, cond, l, delay, nWedge)
%Plots the decay of fluorescent activity when the fly stops moving and the
%re-emergence of activity when the fly starts moving again. Fits an
%exponential to the decay and sigmoidal to the re-emergence.


%%

%l is length of stretch that must be zero
%delay is length of stretch to consider for re-emergence of fluorescence

%find periods of 5+ time points of low vRot+vF
%align beginning and look at fluorescence. Normalize to beginning
%align end and look at fluorescence. normalize to end

threshR = 0;  %want the fly to actually stand still
threshF = 0; %1/100;
pre1 = 2; %dictates how many prior points we include.
pre2 = 1;

name = cond{1}.name;

pers = { 'dark', 'OL', 'CL' };

conds = {'RT', '30C'};

taus = {containers.Map(pers, {[], [], []}) containers.Map(pers, {[], [], []})}; %{RT 30C}

for condid = 1:2

    for per = 1:3

        [DFs, data, flyinds] = getDataShi(cond, pers{per}, 7, 0, nWedge);
        %smooth fluorescence but not vel. Makes it somewhat less noisy
        %nWedges specifies the number of wedges we average over

        maxs = max(DFs{condid}, [], 1);
        
        tPts = data{6*condid-5};
        vR = abs( data{6*condid-4} );
        vF = abs( data{6*condid-3} );
        int = data{6*condid-2};
        mags = data{6*condid-1};
        dirs = data{6*condid};
        vDF = abs(getRates(dirs, tPts, 2*pi) ); %get rates
        flyinds = flyinds{condid}

        

        decays = [];
        delays = [];
        vdecs = [];
        vdels = [];
        vFdecs = [];
        vFdels = [];%magnitude of rate of change of direction of PVA
        magdecs = [];%magnitude of PVA
        magdels = [];
        
        prevind = 1;
        
        xvals1 = [0:l+pre1-1].*(tPts(6)-tPts(5)); %length of decays

        for i = 3:length(vR)-l
            
            %if tPts(i) < tPts(i-1) || i == length(vR)-l %new trial; add taus from previous trial
            if any(flyinds == i) || i == length(vR)-l
                display('starting new fly')
                s = size(decays);
                if s(2) >= prevind %have added bouts of standing
                    mdec = mean(decays(:, prevind:s(2)), 2); %mean intensity

                    fdec = fit(xvals1', mdec, 'exp1');

                    taus{condid}(pers{per}) = [taus{condid}(pers{per}); -fdec.b]; %rate constant

                    prevind = s(2)+1;
                end
            end
            
            %need continuous stationary stretch of length l starting at i
            if all(vR(i:i+(l-1)) <= threshR) && all(vF(i:i+(l-1)) <= threshF ) && ...
                    all(vR(i-pre1:i-1) > threshR) && all(vF(i-pre1:i-1) > threshF ) && ...
                    tPts(i+(l-1)) > tPts(i-pre1)
                
                decays = [decays transpose( int(i-pre1:i+(l-1))/int(i-pre1) ) ]; %normalize to first data point
                
                %{
                if method == 'mean'
                    decays = [decays transpose( int(i-pre1:i+(l-1))/int(i-pre1) ) ]; %normalize to first data point
                elseif method =='max'
                    decays = [decays transpose( maxs(i-pre1:i+(l-1))/maxs(i-pre1) ) ]; %normalize to first data point
                else
                    display('method should be max or mean')
                end
                %}
                
                vFdecs = [vFdecs vDF(i-pre1:i+(l-1)) ]; %no normaliztion
                magdecs = [magdecs transpose( mags(i-pre1:i+(l-1)) ) ];%/ mags(i-pre1) ]; %normalize
                vdecs = [vdecs transpose( vR(i-pre1:i+(l-1)) ) ];
                
                j = l;
                while length(vR) >= i+j && vR(i+j) <= threshR && vF(i+j) <= threshF %find the end of the stationary stretch
                    j = j+1;
                end

                if length(vR) > i+j+(delay-1) %don't want to reach end of list
                    %display('adding delay')
                    %vR( i+j : i+j+(delay-1) )
                    %vF( i+j : i+j+(delay-1) )
                    
                    %want continuous activity
                    if all( vR( i+j : i+j+(delay-1) ) > threshR ) && all( vF( i+j : i+j+(delay-1) ) > threshF ) && ...
                            tPts(i+j+(delay-1)) > tPts(i)
                        %display('still adding delay')
                        
                        delays = [delays transpose( int( i+j-pre2 : i+j+(delay-1) ) / int(i+j-pre2) )];
                        %{
                        if method == 'mean'
                            delays = [delays transpose( int( i+j-pre2 : i+j+(delay-1) ) / int(i+j-pre2) )];
                        elseif method == 'max'
                            delays = [delays transpose( maxs( i+j-pre2 : i+j+(delay-1) ) / maxs(i+j-pre2) )];
                        else
                            display('method should be max or mean')
                        end
                        %}
                        
                        vFdels = [vFdels vDF( i+j-pre2 : i+j+(delay-1) ) ];
                        magdels = [magdels transpose( mags( i+j-pre2 : i+j+(delay-1) ))];% / mags(i+j-pre2) )];
                        vdels = [vdels transpose( vR( i+j-pre2 : i+j+(delay-1) ) ) ];
                    end
                    
                end                
            end
        end

        %vdecs
        %decays
        %vFdecs
        %magdecs
        
        %vdels
        %delays
        %vFdels
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
        
        xvals1 = [0:l+pre1-1].*(tPts(6)-tPts(5)); %length of decays
        xvals2 = [0:delay+pre2-1].*(tPts(6)-tPts(5)); %length of delays
        
        mvdec = mean(vdecs, 2); %mean velocity
        mdec = mean(decays,2); %mean intensity
        mvFdec = mean(vFdecs, 2); %mean rate of change of fluorescence direction
        mmagdec = mean(magdecs, 2); %mean of |PVA|
        
        mvdel = mean(vdels, 2);
        mdel = mean(delays,2);
        mvFdel = mean(vFdels, 2);
        mmagdel = mean(magdels, 2);
        
        newfig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
        
        fdec = fit(xvals1', mdec, 'exp1');
        [fdel, gofdel] = fit(xvals2', mdel, ft, 'problem', mdel(1));     
        fitdec = fdec.a * exp( fdec.b * xvals1); %we fit to exp(kx) NOT exp(-kx)
        fitdel = mdel(1) + fdel.a * ( 1 + exp( -fdel.b * (xvals2 - fdel.c) ) ).^(-1); %fit sigmoidal
        
        fVdec = fit(xvals1', mvFdec, 'exp1');
        [fVdel, gofVdel] = fit(xvals2', mvFdel, ft, 'problem', mvFdel(1));
        fitVdec = fVdec.a * exp( fVdec.b * xvals1 );
        fitVdel = mvFdel(1) + fVdel.a * ( 1 + exp( -fVdel.b * (xvals2 - fVdel.c) ) ).^(-1);
        
        
        fmagdec = fit(xvals1', mmagdec, 'exp1');
        [fmagdel, gofmagdel] = fit(xvals2', mmagdel, ft, 'problem', mmagdel(1));
        fitmagdec = fmagdec.a * exp( fmagdec.b * xvals1 );
        fitmagdel = mmagdel(1) + fmagdel.a * ( 1 + exp( -fmagdel.b * (xvals2 - fmagdel.c) ) ).^(-1);
        
        %plot intensities decay
        subplot(3,2,1)
        hold on
        
        p1 = plot(xvals1, mdec, 'k', 'linewidth', 5)
        title( strcat('intensity-decays','-', name,'-',pers(per),'-',conds{condid}) )
        p2 = plot(xvals1, fitdec, 'b--', 'linewidth', 3)
        legend( { 'Mean', strcat('tau =',num2str(-1/fdec.b, 3)) } )
        plot(xvals1, decays)
        uistack([p2 p1], 'top')
        xlim( [0 max(xvals1)] )

        %plot intensities delay
        subplot(3,2,2)
        hold on
        title( strcat('intensity-delays','-', name, '-',pers(per),'-',conds{condid}) )
        p1 = plot(xvals2, mdel, 'k', 'linewidth', 5)
        p2 = plot(xvals2, fitdel, 'b--', 'linewidth', 3)
        legend( { 'Mean', strcat('delay =',num2str(fdel.c, 3),' slope=',num2str(0.25*fdel.a * fdel.b, 3))    })%,...
            %[num2str(fdel.a) '-' num2str(fdel.b) '-' num2str(fdel.c)] } )
        plot(xvals2, delays)
        uistack([p2 p1], 'top')
        xlim( [0 max(xvals2)] )

        %plot direction decay
        subplot(3,2,3)
        hold on
        p1 = plot(xvals1, mvFdec, 'k', 'linewidth', 5)
        p2 = plot(xvals1, mvdec, 'k--', 'linewidth', 3)
        title( strcat('rate-decays','-', name,'-',pers(per),'-',conds{condid}) )
        p3 = plot(xvals1, fitVdec, 'b--', 'linewidth', 3)
        legend( { 'Mean', 'vRot', strcat('tau =',num2str(-1/fVdec.b, 3)) } )
        plot(xvals1, vFdecs)
        uistack([p3 p2 p1], 'top')
        xlim( [0 max(xvals1)] )

        %plot direction delay
        subplot(3,2,4)
        hold on
        p1 = plot(xvals2, mvFdel, 'k', 'linewidth', 5)
        p2 = plot(xvals2, mvdel, 'k--', 'linewidth', 3)
        p3 = plot(xvals2, fitVdel, 'b--', 'linewidth', 3)
        title( strcat('rate-delays','-', name, '-',pers(per),'-',conds{condid}) )
        legend( { 'Mean', 'vRot', strcat('delay =',num2str(fVdel.c, 3),' slope=',num2str(0.25*fVdel.a * fVdel.b, 3))    })%,...
            %[num2str(fdel.a) '-' num2str(fdel.b) '-' num2str(fdel.c)] } )
        plot(xvals2, vFdels)
        uistack([p3 p2 p1], 'top')
        xlim( [0 max(xvals2)] )
        
        
        %plot |PVA| decay
        subplot(3,2,5)
        hold on
        p1 = plot(xvals1, mmagdec, 'k', 'linewidth', 5)
        p2 = plot(xvals1, fitmagdec, 'b--', 'linewidth', 3)
        plot(xvals1, magdecs)
        title( strcat('|PVA|-decays','-', name,'-',pers(per),'-',conds{condid}) )
        legend( { 'Mean', strcat('tau =',num2str(-1/fmagdec.b, 3)) } )
        uistack([p2 p1], 'top')
        xlim( [0 max(xvals1)] )

        %plot |PVA| delay
        subplot(3,2,6)
        hold on
        p1 = plot(xvals2, mmagdel, 'k', 'linewidth', 5)
        p2 = plot(xvals2, fitmagdel, 'b--', 'linewidth', 3)
        plot(xvals2, magdels)
        title( strcat('|PVA|-delays','-', name, '-',pers(per),'-',conds{condid}) )
        legend( { 'Mean', strcat('delay =',num2str(fmagdel.c, 3),' slope=',num2str(0.25*fmagdel.a * fmagdel.b, 3))    })%,...
            %[num2str(fdel.a) '-' num2str(fdel.b) '-' num2str(fdel.c)] } )
        uistack([p2 p1], 'top')
        xlim( [0 max(xvals2)] )
        
        
        
        newfig.PaperUnits = 'inches';
        newfig.PaperPosition = [0 0 8 11.5];
        
        %dest = strcat( dir, 'decay_delay_', name, '_', pers(per), '_', conds{condid} )
        dest = strcat(dir, 'decayDelay_', name, '_', conds{condid}, '_', pers{per}, '_', 'nWedge', num2str(nWedge));
        
        print(newfig, dest, '-dpdf');


    end


end

