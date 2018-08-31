%function headingVectorMag(a,b, period)

%make plot of heading vector magnitude vs. windowlength
%do this for all flies and plot mean =- standard deviation for four trials
%of interest on same plot.
%at each timepoint, calculate percent change and whether or not it's
%siginificant cf. control.

periods = {'dark', 'CL', 'OL'};

data = { [] [] [] [] [] [] [] [] [] [] [] [] };
DFs = {[] []};
conds = {'All' 'All_30C'};
names = {'RT', '30C'};
smoothF = 0;
smoothV = 0;




dirs = { '/Users/loaner/Documents/imaging/Data_Dan/shi/empty/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/EPG/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/GE/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEG/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEN2/'}%,...
    %'/Users/loaner/Documents/imaging/Data_Dan/shi/PEN1/'};

for j = 1:2; %iterate over conditions
        
    means0 = cell(1,3)
    stdevs0 = cell(1,3)
    
    for d = 1:length(dirs)
        dir = dirs{d}
        
        try
            from_file = load(strcat(dir, 'cond'), 'cond');
            cond=from_file.cond;    
        catch
            cond = FlyDatLoad(1);
            save(strcat(dir, 'cond'), 'cond');
        end
        
        name = cond{1}.name;
    
        fig1 = figure('units','normalized','outerposition',[0 0 1 1])%, 'visible', 'off');
        
        fig2 = figure('units','normalized','outerposition',[0 0 1 1])%, 'visible', 'off');
        

        
        for perid = 1:3
            
            period = periods{perid};
    
            mags = [];
            for i = 1:length(cond{1}.allFlyData); %iterate over flies

                for k = 1:2 %iterate over trials



                    if j == 1
                        ind = k;
                    else
                        ind = k+1; %don't use first trial from 30C
                    end

                    [tPts,darkPer, OLPer, CLPer, CWPer, CCWPer, DF, heading, headingPlt, vRot, vF, stripePos, stripePosPlt, stripeJumps]...
                        = extractShiData(cond, i, conds{j}, ind, smoothF, smoothV);

                    %get indices of period of interest
                    if strcmp(period, 'dark')
                        per = darkPer;
                    elseif strcmp(period, 'OL')
                        per = OLPer;
                    elseif strcmp(period, 'CL')
                        per = CLPer;
                    else
                        display('per not recognized, please select dark or OL')
                    end

                    if i == 1 && k ==1 && d==1 && perid == 0 %first trial
                        windows = 1:length(per)-3; %give us some slack for subseqeuent periods
                        L = length(windows);
                        ts = tPts(per);
                        ts = ts(1:L)-ts(1);
                    end

                    heading = heading(per);

                    vectormags = zeros(1, L);

                    for window = 1:L

                        vectormags(window) = meanvecmag(heading, window);
                    end

                    mags = [mags vectormags'];
                end
            end

            means = mean(mags, 2);
            stdevs = std(mags, 0, 2);
            
            if d == 1 %store shibire for future comparison
                means0{perid} = means;
                stdevs0{perid} = stdevs;
            end

            figure(fig1)
            subplot(3,1, perid)
            hold on
            x = [ts',fliplr(ts')];
            yy = [means'-stdevs',fliplr(means'+stdevs')];
            fill(x,yy,'g','facealpha',.1, 'LineStyle', 'none')
            plot(ts, means, 'g')
            %plot(ts, means+stdevs, 'b:')
            %plot(ts, means-stdevs, 'b:')
            ylabel('Mean heading vector magnitude')
            xlabel('Window length (s)')
            ylim([0 1])
            xlim([ts(1) ts(end)])
            if j ==1
                title(strcat(name, '-RT-', period))
            else
                title(strcat(name, '-30C-', period))
            end
            
            figure(fig2)
            subplot(3,1, perid)
            
            hold on
            yy = [means0{perid}'-stdevs0{perid}',fliplr(means0{perid}'+stdevs0{perid}')];
            fill(x,yy,'r','facealpha',.1, 'LineStyle', 'none')
            yy = [means'-stdevs',fliplr(means'+stdevs')];
            fill(x,yy,'g','facealpha',0.07, 'LineStyle', 'none')
            plot(ts, means0{perid}, 'r')
            plot(ts, means, 'g')
            ylabel('Mean heading vector magnitude')
            xlabel('Window length (s)')
            ylim([0 1])
            xlim([ts(1) ts(end)])
            legend({'control', name})
            if j ==1
                title(strcat(name, '-RT-', period))
            else
                title(strcat(name, '-30C-', period))
            end

        end
        figure(fig1)
        fig1.PaperUnits = 'inches';
        fig1.PaperPosition = [0 0 8 11.5];
        %print(fig1, 'test', '-dpdf')
        print(fig1, strcat(dir, 'headingVectorMag/headingVectorMag_', name, '_', names{j}), '-dpdf'); %save
        figure(fig2)
        fig2.PaperUnits = 'inches';
        fig2.PaperPosition = [0 0 8 11.5];
        print(fig2, strcat(dir, 'headingVectorMag/headingVectorMag_compare_', name, '_', names{j}), '-dpdf'); %save
    end
end

                