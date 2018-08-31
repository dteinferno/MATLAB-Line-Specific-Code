function [pres, posts] = persistance(dir, cond, l)
%Consider correlation between PVA before standing still and when
%re-initiating movement
%Do this across all flies because we don't have enough data?
%alternatively on a per-fly basis and consider mean correlation

%%

%l is length of stretch that must be zero to consider it a bout of standing
%still


threshR = 0;  %want the fly to actually stand still
threshF = 0;


name = cond{1}.name;

pers = { 'dark', 'OL', 'CL' };

conds = {'RT', '30C'};

for condid = 1:2
    newfig = figure('units','normalized','outerposition',[0 0 1 1])%, 'visible', 'off');
    for per = 1:3

        [DFs, data] = getDataShi(cond, pers{per}, 0, 0); %can't smoothen for this analysis

        tPts = data{6*condid-5};
        vR = abs( data{6*condid-4} );
        vF = abs( data{6*condid-3} );
        int = data{6*condid-2};
        mags = data{6*condid-1};
        dirs = data{6*condid}-pi;
        vDF = abs(getRates(dirs, tPts, 2*pi) ); %get rates

        vFs = [];
        vRs = [];%magnitude of rate of change of direction of PVA
        %vR2s = [];
        %vR3s = [];
        pres = [];
        posts = [];
        ts = [];

        for i = 3:length(vR)-(l)
            %need continuous stationary stretch of length l starting at i
            if all(vR(i:i+(l-1)) <= threshR) && all(vF(i:i+(l-1)) <= threshF ) && ... %stand still for l frames
                    all(vR(i-1) > threshR) && all(vF(i-1) > threshF &&...
                    tPts(i+(l-1)) > tPts(i-1) ) %need a new standing still bout
                
                j = l;
                
                try
                    while (vR(i+j) <= threshR && vF(i+j) <= threshF) || mags(i+j) <0.2 
                        %find the end of the stationary stretch and require
                        %|PVA| > 0.2
                        j = j+1;
                    end
                    
                    if tPts(i+j) > tPts(i) %don't want to start new trial
                    
                        pres = [pres dirs(i-1)]; %last moving fame
                        posts = [posts dirs(i+j)]; %first moving frame
                        vRs = [vRs [vR(i-1);vR(i);vR(i+j)] ];

                        %vFs = [vFs transpose( vF(i-1:i+j+1) ) ];
                        ts = [ts tPts(i+j)-tPts(i)];
                    end
                catch %we run out of list
                end
                    
            end
        end

        vRs
        pres
        posts
        
        difs = getOffset(posts, pres, 2*pi); %get closest distance to posts from pres
        posts_reg = pres+difs; %this allows us to 'wrap around' for the regression

        
        subplot(3,2,2*per-1)
        
        scatter(pres, posts_reg)
        hold on
        
        p1 = polyfit(pres,posts_reg,1);
        yfit1 = polyval(p1,pres);
        yresid = posts_reg - yfit1;
        SSresid = sum(yresid.^2);   
        SStotal = (length(posts_reg)-1) * var(posts_reg);
        rsq1 = 1 - SSresid/SStotal;
        
        rho = corr(transpose(pres), transpose(posts_reg)); %measure correlation
        
        plot(pres, yfit1, 'r')
        
        title( sprintf(strcat(name, '-', conds{condid}, '-', pers{per}, ' slope=%.2f R=%.2f'), p1(1), rho ) )
        ylabel('post')
        xlabel('pre')
        xlim([-pi, pi])
        ylim([1.05*min(posts_reg), 1.05*max(posts_reg)])
        
        
        
        subplot(3,2,2*per)
        
        scatter(ts, difs)
        hold on
        
        p1 = polyfit(ts,difs,1);
        yfit1 = polyval(p1,ts);
        yresid = difs - yfit1;
        SSresid = sum(yresid.^2);   
        SStotal = (length(difs)-1) * var(difs);
        rsq1 = 1 - SSresid/SStotal;
        
        rho = corr(transpose(ts), transpose(difs)); %measure correlation
        
        plot(ts, yfit1, 'r')
        
        title( sprintf(strcat(name, '-', conds{condid}, '-', pers{per}, ' slope=%.2f R=%.2f'), p1(1), rho ) )

        xlabel('time')
        ylabel('PVA post - pre')
        
        xlim( [0 1.05*max(ts) ] )
        
        
        



    end
    
    newfig.PaperUnits = 'inches';
    newfig.PaperPosition = [0 0 8 11.5];
    print(newfig, strcat(dir, strcat('persistance/persistance_', name, '_', conds{condid})), '-dpdf'); %save violinplot


end

