function data = scatterVelEB1C(dir, cond, period, nWedge)

%Makes scatterplots and violin plots for intensity and PVA magnitude against
%rotational velocity for dark, CL or OL visual stimulus specified by
%'period'. Saves to dir. Pools over flies and trials.

vels = [0 0 pi/6 pi/3 2*pi/3 20]; %rotational velocity bins

name = cond{1}.name;

data = { [] [] [] [] [] [] };
names = {'vRot' strcat('intensity ', name, ' RT') strcat('|PVA| ', name, ' RT')... 
'vRot' strcat('intensity ', name, ' 30C') strcat('|PVA| ', name, ' 30C')};
    
conds = {'All' 'All_30C'};

%% Get Data

for i = 1:length(cond{1}.allFlyData); %iterate over flies
    
    for j = 1:2; %iterate over conditions (hot and cold)
   
        for k = 1:2 %iterate over trials
            
            if j == 1
                ind = k;
            else
                ind = k+1; %don't use first trial from 30C
            end
            
            %smooth with SG filter since we will be looking at directions
            [tPts,darkPer, OLPer, CLPer, CWPer, CCWPer, DF, heading, headingPlt, vRot, vF, stripePos, stripePosPlt, stripeJumps]...
                = extractShiData(cond, i, conds{j}, ind, 5, 5); 
            
            int = mean( maxk(DF, nWedge, 1), 1)
            
            %getmagnitude of PVAs
            m = zeros(1, length(int));
            for ind = 1:length(int);
                vec = DF(:, ind);
                [~, mi] = getVecSum( vec );
                m(ind) = mi/sum(abs(vec)); %normalize |PVA|
            end
            
            vRot = transpose(abs(vRot)); %work with row vectors
            
            %get period of interest
            if strcmp(period, 'dark')
                per = darkPer;
            elseif strcmp(period, 'OL')
                per = OLPer;
            elseif strcmp(period, 'CL')
                per = CLPer;
            else
                display('per not recognized, please select dark or OL')
            end
            
            data{3*j-2} = [data{3*j-2} vRot(per)]; %1-3 is cold, 4-6 hot
            data{3*j-1} = [data{3*j-1} int(per)]; %intensity
            data{3*j} = [data{3*j} m(per)]; %|PVA|
            
        end
    end
end
    
%% Plot scatterplots

fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

%iterate over hot and cold
for k=1:2
    subplot(2, 2, 2*k-1)

    scatter(data{3*k-2}, data{3*k-1}, 3, 'filled') %vRot vs. itensity

    xlabel(names{3*k-2})
    ylabel(names{3*k-1})


    subplot(2, 2, 2*k)

    scatter(data{3*k-2}, data{3*k}, 3, 'filled') %vRot vs |PVA|

    xlabel(names{3*k-2})
    ylabel(names{3*k})

end         
print(fig, strcat(dir, 'scatter_intensity', name, '_', period, '_nWedge', num2str(nWedge)), '-dpdf');

%%  Plot Boxplots

box = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off'); %plot violin plots

set(box, 'DefaultTextFontSize', 4);


for k=1:2 %iterate over hot and cold
    
    for j = 1:2 %iterate over m, |PVA|
        subplot(4, 1, 2*k-2+j)
        
        if j == 1
            ymax = 5.0;
        else
            ymax = 1.0;
        end

        boxdat = [];
        labels = {};
        xlabs = [];
        distDat = [];%
        for m = 1:length(vels)-1
            if m == 1
                dat = data{3*k-2+j}( data{3*k-2} == 0 ); %zero velocity
                dat0 = dat; %reference data at zero vel for stats and normalization
                ps = 1;
                groups = zeros(length(dat), 1); %specify groupings
            else
                dat = data{3*k-2+j}( data{3*k-2} > vels(m) & data{3*k-2} <= vels(m+1) ); %get data for bin of interest
                if length(dat) > 0
                    groups = vertcat( groups, (m-1)*ones(length(dat),1) ); %add to groupings if we have data
                end

            end

            if length(dat) > 0
                vertcat(distDat, dat);%
                [~, p] = ttest2(dat, dat0, 'Tail', 'both', 'Vartype', 'unequal'); %check vs ref data
                xlabs = [xlabs length(dat) p]; %add some data to xlabels
                if j == 1
                    dat = dat / mean(dat0); %normalize intensity to zero velocity
                end
                boxdat = [boxdat dat];
                M = mean(dat);
                dev = std(dat);
                labels{m} = sprintf('%.2f m=%.2f(%.2f)', vels(m+1), M, dev); %add some data to labels
            end
        end

        xlab = '';
        for m = 1:length(ps)
            xlab = strcat(xlab, ',  n=%d p=%.2e'); %add some more data to labels
        end

        distributionPlot(transpose(boxdat), 'groups', transpose(groups), 'color', 'b',...
            'showMM', 5, 'xNames', labels) %requires distributionPlot package
        ylim( [0 ymax] )
        set(gca,'FontSize',8);
        %boxplot(boxdat, groups, 'Labels',labels, 'Whisker', 5)
        %set(findobj(gca,'Type','text'),'FontSize',18)
        xlabel( sprintf(xlab, xlabs), 'fontsize', 9)
        title(names{3*k-2+j}, 'FontSize', 16)
    end
end
box.PaperUnits = 'inches';
box.PaperPosition = [0 0 8 11.5];
print(box, strcat(dir, 'box_intensity', name, '_', period, '_nWedge', num2str(nWedge)), '-dpdf');
%}