function scatterVelEB(dir, green, red)
%given a directory in which to dump files and names of green & red channels
%(i.e. celltype),  makes scatterplots and violinplots of closedloop data
%(vR, vF, intensity, |PVA|). Requires the 'distributionPlot' package.

try %if we've already run FlyDatLoad we load the previous data
    from_file = load(strcat(dir, 'cont'), 'alldata');
    alldata=from_file.alldata;
    
catch %if not, we do it now and store data
    alldata = FlyDatLoad(2, 'EB');
    save(strcat(dir, 'cont'), 'alldata');
end



cd ~/Documents/Imaging/Data_Dan

vels = { [0 0 pi/6 pi/3 2*pi/3 20], [0 0 0.2 0.5 1 20]}; %bin rotational and forward velocities

for i = 1:length(alldata{1}.allFlyData); %iterate over flies

    fly = alldata{1}.allFlyData{i};
    %we pool data by fly
    data = { [] [] [] [] [] [] };
    names = {'vRot' 'vF' strcat('intensity ', green) strcat('intensity ', red),...
    strcat('|PVA| ', green) strcat('|PVA| ', red) }; %for later use in labelling diagrams
    
    try%we may not have consistent nomenclature
        L = length(fly.Dark);
    catch
        L = length(fly.All);
    end
    
    for j = 1:L; %iterate over trials

        try
            trial = fly.Dark{j}
        catch
            trial = fly.All{j}
        end
   
        if length(trial) > 0 & max(trial.positionDatMatch.vF) > 0 %discard trials where fly stands still throughout
            
            %get imaging and velocity data
            datG = trial.GROIaveMax-1; %numbering is counterclockwise?
            datR = trial.RROIaveMax-1; 
            vR = trial.positionDatMatch.vRot( trial.positionDatMatch.Closed(1:length(trial.positionDatMatch.vRot))== 1 );
            %discard open loop data for now
            vF = trial.positionDatMatch.vF( trial.positionDatMatch.Closed(1:length(trial.positionDatMatch.vRot))== 1 );
            
            smooth = 3 %smoothing by averaging at this stage, have gone over to SG smoothing, may want to rewrite
            if smooth > 0
                s = size(datG);
                s = s(1);
                display('smoothening data');
                vR = transpose(Smooth(vR, smooth)); %work with row vectors
                vF = transpose(Smooth(vF, smooth));
                newG = [];
                newR = [];
                for ind = 1:s;
                    g = datG(ind,:);
                    newG(ind,:) = Smooth(g(trial.positionDatMatch.Closed== 1), smooth);
                    r = datR(ind,:);
                    newR(ind,:) = Smooth(r(trial.positionDatMatch.Closed== 1), smooth);
                end
                datG = newG;
                datR = newR;    
            end
            
            intG = mean(datG,1); %get the mean intensities (mean dF/F across ROIs)
            intR = mean(datR,1);
            
            %getmagnitude of PVAs
            mG = zeros(1, length(intG));
            mR = zeros(1, length(intG));
            for ind = 1:length(intG);
                vecG = datG(:, ind);
                vecR = datR(:, ind);
                
                [~, mGi] = getVecSum( vecG );
                mG(i) = mGi/sum(abs(vecG)); %add normalized PVA magnitude

                [~, mRi] = getVecSum( vecR );
                mR(i) = mRi/sum(abs(vecR));
            end
            
            vR = abs(vR); %looking at EB, use magnitude of velocities
            vF = abs(vF);
            intG = intG(1:length(vR)); %only take data for which we have velocities
            intR = intR(1:length(vR));
            mG = mG(1:length(vR));
            mR = mR(1:length(vR));
            
            %add data for trial to pooled data for fly
            data{1} = [data{1} vR];
            data{2} = [data{2} vF];
            data{3} = [data{3} intG];
            data{4} = [data{4} intR];
            data{5} = [data{5} mG];
            data{6} = [data{6} mR];
            
        end
    end
    
    %% Plot scatterplots
    name = strcat(dir, sprintf('fly%d_', i)); %save to dict labelled by flyID
    for j=1:2 %iterate over vRot and vF on x axis
        
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off') %make scatterplot

        %iterate over green and red
        for k=1:2
            subplot(2, 2, 2*k-1)
            
            scatter(data{j}, data{2*k+1}, 3, 'filled') %intensities
            
            xlabel(names{j})
            ylabel(names{2*k+1})
            
            
            subplot(2, 2, 2*k)
            
            scatter(data{j}, data{2*k+2}, 3, 'filled') %magnitudes
            
            xlabel(names{j})
            ylabel(names{2*k+2})
            
        end         
        print(fig, strcat(name, strcat('scatter_', names{j})), '-dpdf'); %save to dict
        
        box = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off') %make boxplot

        for k=1:4 %iterate over (vRot, vF) x (grenn, red)
            subplot(4, 1, k)
            
            boxdat = [];
            labels = {};
            xlabs = [];
            distDat = [];%
            for m = 1:length(vels{j})-1
                if m == 1
                    dat = data{k+2}( data{j} == 0 ); %find staitonary data
                    dat0 = dat;% store stationary data for comparison
                    ps = [1];
                    groups = zeros(length(dat), 1); %group for plotting
                else
                    dat = data{k+2}( data{j} > vels{j}(m) & data{j} <= vels{j}(m+1) ); %get data in velocity bin
                    if length(dat) > 0
                        groups = vertcat( groups, (m-1)*ones(length(dat),1) ); %add to groupings
                    end
                    
                end
                
                if length(dat) > 0
                    vertcat(distDat, dat);%
                    [r, p] = ttest2(dat, dat0, 'Tail', 'both', 'Vartype', 'unequal') %run ttest relative to stationary data
                    xlabs = [xlabs length(dat) p]; %construct xlabels
                    if k==1 || k==2
                        dat = dat / mean(dat0); %normalize intensities
                    end
                    boxdat = [boxdat dat];
                    M = mean(dat);
                    dev = std(dat);
                    labels{m} = sprintf('%.2f m=%.2f(%.2f)', vels{j}(m+1), M, dev); %update labels with some info on the data
                end
            end

            xlab = ''
            for m = 1:length(ps)
                xlab = strcat(xlab, ',  n=%d p=%.2e') 
            end

            distributionPlot(transpose(boxdat), 'groups', transpose(groups), 'color', 'b',...
                'showMM', 5, 'xNames', labels) %make violin plot
            ylim( [0 max(boxdat)] )
            set(gca,'FontSize',8);
            xlabel( sprintf(xlab, xlabs), 'fontsize', 9)
            title(names{k+2}, 'FontSize', 16) 
            
        end
        box.PaperUnits = 'inches';
        box.PaperPosition = [0 0 8 11.5];
        print(box, strcat(name, strcat('violin_', names{j})), '-dpdf'); %save violinplot
    end
end