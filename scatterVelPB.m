function scatterVelPB(dir, green, red)
%given a directory in which to dump files and names of green & red channels (i.e. celltype), 
%makes scatterplots and violinplots of closedloop data (vR, vF, intensity, |PVA| for left and right PB).
%Requires the 'distributionPlot' package. If FlyDatLoad has already been run and the container stored in dir, we load that
% otherwise, we rerun FlyDatLoad. Was written for FlyDatLoad_old, have not tested with new FlyDatLoad.

try %load data from dir if possible
    from_file = load(strcat(dir, 'cont'), 'alldata');
    alldata=from_file.alldata;    
catch
    alldata = FlyDatLoad(2);
    save(strcat(dir, 'cont'), 'alldata');
end

vels = [-20 -pi/2 -pi/6 0 0 pi/6 pi/2 20]; %bin rotational velocities
dirs = {'L' 'R'};
names = {'vR' { strcat('intensity ', green) strcat('intensity ', red) ...
    strcat('|PVA|', green) strcat('|PVA|', red) } ...
    { strcat('intensity ', green) strcat('intensity ', red) ...
    strcat('|PVA|', green) strcat('|PVA|', red) } }; %for later use in labelling diagrams
    
for i = 1:length(alldata{1}.allFlyData);

    fly = alldata{1}.allFlyData{i};
            
    data = { [] {[] [] [] []} {[] [] [] []} }; %pool data by fly
    
    try
        L = length(fly.Dark);
    catch
        L = length(fly.All);
    end
    
    for j = 1:L; %iterate over trials

        try
            trial = fly.Dark{j};
        catch
            trial = fly.All{j};
        end
   
        if length(trial) > 0 & max(trial.positionDatMatch.vF) > 0
            
            datG = trial.GROIaveMax-1; %CCW numbering?
            datR = trial.RROIaveMax-1; 
            %consider only closed loop data
            vR = trial.positionDatMatch.vRot( trial.positionDatMatch.Closed(1:length(trial.positionDatMatch.vRot))== 1 );
            
            %Should change to SG filtering
            smooth = 3;
            if smooth > 0
                s = size(datG);
                s = s(1);
                display('smoothening data');
                vR = transpose(Smooth(vR, smooth)); %work with row vectors
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
            
            for n = 1:2 %iterate over left and right PB
                
                intG = sum(datG(n*9-8:n*9,:),1);
                intR = sum(datR(n*9-8:n*9,:),1);
                
            %getmagnitude of PVAs
                mG = zeros(1, length(intG));
                mR = zeros(1, length(intG));
                for ind = 1:length(intG);
                    vecG = datG(n*9-8:n*9, ind);
                    vecR = datR(n*9-8:n*9, ind);

                    [~, mGi] = getVecSum( vecG );
                    mG(ind) = mGi/sum(abs(vecG)); %add PVA magnitude (normalized)

                    [~, mRi] = getVecSum( vecR );
                    mR(ind) = mRi/sum(abs(vecR));
                end
                intG = intG(1:length(vR)); %only consider data for which we have vRot
                intR = intR(1:length(vR));
                mG = mG(1:length(vR));
                mR = mR(1:length(vR));

                data{n+1}{1} = [data{n+1}{1} intG]; %update our cumulative data
                data{n+1}{2} = [data{n+1}{2} intR];
                data{n+1}{3} = [data{n+1}{3} mG];
                data{n+1}{4} = [data{n+1}{4} mR];
            end

            data{1} = [data{1} vR];
            
        end
        
    end
    
    %% Plot scatterplots
    
    name = strcat(dir, sprintf('fly%d_', i));
    for j=1:2 %iterate over left, right
 
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
        
        %iterate over green and red
        for k=1:2
            subplot(2, 2, 2*k-1)
            
            scatter(data{1}, data{j+1}{2*k-1}, 5, 'filled')
            
            xlabel(names{1})
            ylabel(names{j+1}{2*k-1})
            
            
            subplot(2, 2, 2*k)
            
            scatter(data{1}, data{j+1}{2*k}, 5, 'filled')
            
            xlabel(names{1})
            ylabel(names{j+1}{2*k})
            
        end         
        print(fig, strcat(name, strcat('scatter_', dirs{j})), '-dpdf');
        
        %% Make violin plots / boxplots

        box = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
        
        for k=1:4 %iterate over red and green, intensity, |PVA|
            subplot(4, 1, k);
            
            dat0 = data{j+1}{k}( data{1} == 0 ); %standing still
            boxdat = [];
            labels = {};
            xlabs = [];
            groups = [];
            n = 1;
            for m = 1:length(vels)-1

                if vels(m) == 0 && vels(m+1) == 0
                    dat = dat0; %standing stell
                else
                    dat = data{j+1}{k}( data{1} > vels(m) & data{1} <= vels(m+1) ); %bin by velocity
                end
                if length(dat) > 0
                    groups = vertcat( groups, (m-1)*ones(length(dat),1) ); %add to groups

                    [r, p] = ttest2(dat, dat0, 'Tail', 'both', 'Vartype', 'unequal') %run stats relative to standing still
                    xlabs = [xlabs length(dat) p];
                    if k==1 || k==2 %normalize intensities
                        dat = dat / mean(dat0);
                    end
                    boxdat = [boxdat dat];
                    M = mean(dat);
                    dev = std(dat);
                    labels{n} = sprintf('%.2f m=%.2f(%.2f)', vels(m+1), M, dev); %put some data in labels
                    n = n+1; %number the non-empty groups
                end
            end

            xlab = ''
            for m = 1:length(labels)
                xlab = strcat(xlab, ', n=%d p=%.2e'); %add some labels
            end

            distributionPlot(transpose(boxdat), 'groups', transpose(groups), 'color', 'b',...
                'showMM', 5, 'xNames', labels) %plot violinplot
            ylim( [0 max(boxdat)] )
            %boxplot(boxdat, groups, 'Labels',labels, 'Whisker', 5)
            set(gca,'FontSize',7);
            xlabel( sprintf(xlab, xlabs), 'fontsize', 8)
            title(names{j+1}{k}, 'FontSize', 16) 
            
        end
        box.PaperUnits = 'inches';
        box.PaperPosition = [0 0 8 11.5];
        %set(findobj(gca,'Type','text'),'fontsize',9);
        print(box, strcat(name, strcat('violin_', dirs{j})), '-dpdf'); %save
    end
end