function scatterVelPB(dir, green, red)

cd ~/Documents/Imaging/Data_Dan

%dir = '~/Documents/Imaging/Data_Dan/PEN1_G_EPG_R_PB/';
%green = 'PEN1';
%red = 'EPG';

%dir = '~/Documents/Imaging/Data_Dan/PEN2_G_PEG_R_PB/';
%green = 'PEN2';
%red = 'PEG';

%dir = '~/Documents/Imaging/Data_Dan/PEN2_R_EPG_G_PB/';
%green = 'EPG';
%red = 'PEN2';

%dir = '~/Documents/Imaging/Data_Dan/D7_R_EPG_G_PB/';
%green = 'EPG';
%red = 'D7';

%dir = '~/Documents/Imaging/Data_Dan/D7_G_EPG_R_PB/';
%green = 'D7';
%red = 'EPG';

try
    from_file = load(strcat(dir, 'cont'), 'alldata');
    alldata=from_file.alldata;    
catch
    alldata = FlyDatLoad(2, 'EB');
    save(strcat(dir, 'cont'), 'alldata');
end

vels = [-20 -pi/2 -pi/6 0 0 pi/6 pi/2 20];
dirs = {'L' 'R'};
names = {'vR' { strcat('intensity ', green) strcat('intensity ', red) ...
    strcat('|PVA|', green) strcat('|PVA|', red) } ...
    { strcat('intensity ', green) strcat('intensity ', red) ...
    strcat('|PVA|', green) strcat('|PVA|', red) } };
    
for i = 1:length(alldata{1}.allFlyData);

    fly = alldata{1}.allFlyData{i};
            
    data = { [] {[] [] [] []} {[] [] [] []} };
    
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
            
            datG = trial.GROIaveMax-1; %which way around does the numbering go???
            datR = trial.RROIaveMax-1; 
            vR = trial.positionDatMatch.vRot( trial.positionDatMatch.Closed(1:length(trial.positionDatMatch.vRot))== 1 );
            
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
            
            for n = 1:2
                
                intG = sum(datG(n*9-8:n*9,:),1);
                intR = sum(datR(n*9-8:n*9,:),1);
                
            %getmagnitude of PVAs
                mG = []
                mR = []
                for ind = 1:length(intG);
                    vecG = datG(n*9-8:n*9, ind);
                    vecR = datR(n*9-8:n*9, ind);

                    [dGi, mGi] = getVecSum( vecG );
                    mG = [mG, mGi/sum(abs(vecG))];

                    [dRi, mRi] = getVecSum( vecR );
                    mR = [mR, mRi/sum(abs(vecR))];
                end
                intG = intG(1:length(vR))
                intR = intR(1:length(vR))
                mG = mG(1:length(vR))
                mR = mR(1:length(vR))

                data{n+1}{1} = [data{n+1}{1} intG]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
                data{n+1}{2} = [data{n+1}{2} intR]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
                data{n+1}{3} = [data{n+1}{3} mG]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
                data{n+1}{4} = [data{n+1}{4} mR]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
            end

            data{1} = [data{1} vR]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
            
        end
        
    end
    
    %% Plot scatterplots
    
    name = strcat(dir, sprintf('fly%d_', i));
    for j=1:2 %iterate over left, right
        if i == 1
            fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
        else
            fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
        end
        
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
        
        %% Make violing plots / boxplots
        if i == 1
            box = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
        else
            box = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
        end
        %set(box, 'DefaultTextFontSize', 4);
        
        for k=1:4 %iterate over red and green
            subplot(4, 1, k)
            
            dat0 = data{j+1}{k}( data{1} == 0 )        
            boxdat = [];
            labels = {};
            xlabs = [];
            groups = []
            n = 1
            for m = 1:length(vels)-1

                if vels(m) == 0 & vels(m+1) == 0
                    dat = dat0
                else
                    dat = data{j+1}{k}( data{1} > vels(m) & data{1} <= vels(m+1) );
                end
                    if length(dat) > 0
                    groups = vertcat( groups, (m-1)*ones(length(dat),1) );
                end
                
                if length(dat) > 0
                    [r, p] = ttest2(dat, dat0, 'Tail', 'both', 'Vartype', 'unequal')
                    xlabs = [xlabs length(dat) p];
                    if k==1 || k==2
                        dat = dat / mean(dat0);
                    end
                    boxdat = [boxdat dat];
                    M = mean(dat);
                    dev = std(dat);
                    labels{n} = sprintf('%.2f m=%.2f(%.2f)', vels(m+1), M, dev);
                    n = n+1
                end
            end

            xlab = ''
            for m = 1:length(labels)
                xlab = strcat(xlab, ', n=%d p=%.2e') 
            end
            
            newlabels = {}
            distributionPlot(transpose(boxdat), 'groups', transpose(groups), 'color', 'b',...
                'showMM', 5, 'xNames', labels)
            ylim( [0 max(boxdat)] )
            %boxplot(boxdat, groups, 'Labels',labels, 'Whisker', 5)
            set(gca,'FontSize',7);
            xlabel( sprintf(xlab, xlabs), 'fontsize', 8)
            title(names{j+1}{k}, 'FontSize', 16) 
            
        end
        box.PaperUnits = 'inches';
        box.PaperPosition = [0 0 8 11.5];
        %set(findobj(gca,'Type','text'),'fontsize',9);
        print(box, strcat(name, strcat('violin_', dirs{j})), '-dpdf');
    end
end