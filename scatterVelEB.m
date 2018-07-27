function scatterVelEB(dir, green, red)

%dir = '~/Documents/Imaging/Data_Dan/PEN2_G_PEG_R_EB/';
%green = 'PEN2';
%red = 'PEG';

%dir = '~/Documents/Imaging/Data_Dan/PEN1_G_EPG_R_EB/';
%green = 'PEN1';
%red = 'EPG';

%dir = '~/Documents/Imaging/Data_Dan/EPG_G_GE_R_EB/';
%green = 'EPG';
%red = 'GE';

%dir = '~/Documents/Imaging/Data_Dan/PEN2_G_EPG_R_EB/';
%green = 'PEN2';
%red = 'EPG';

%dir = '~/Documents/Imaging/Data_Dan/PEN2_R_EPG_G_EB/';
%green = 'EPG';
%red = 'PEN2';

try
    from_file = load(strcat(dir, 'cont'), 'alldata');
    alldata=from_file.alldata;
    
catch

    alldata = FlyDatLoad(2, 'EB');
    save(strcat(dir, 'cont'), 'alldata');
end



cd ~/Documents/Imaging/Data_Dan

vels = { [0 0 pi/6 pi/3 2*pi/3 20], [0 0 0.2 0.5 1 20]};

for i = 1:length(alldata{1}.allFlyData);

    fly = alldata{1}.allFlyData{i};
            
    data = { [] [] [] [] [] [] };
    names = {'vRot' 'vF' strcat('intensity ', green) strcat('intensity ', red),...
    strcat('|PVA| ', green) strcat('|PVA| ', red) };
    
    try
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
   
        if length(trial) > 0 & max(trial.positionDatMatch.vF) > 0
            
            datG = trial.GROIaveMax-1; %which way around does the numbering go???
            datR = trial.RROIaveMax-1; 
            vR = trial.positionDatMatch.vRot( trial.positionDatMatch.Closed(1:length(trial.positionDatMatch.vRot))== 1 )
            vF = trial.positionDatMatch.vF( trial.positionDatMatch.Closed(1:length(trial.positionDatMatch.vRot))== 1 )
            
            smooth = 3
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
            
            intG = sum(datG,1)
            intR = sum(datR,1)
            
            %getmagnitude of PVAs
            mG = []
            mR = []
            for ind = 1:length(intG);
                vecG = datG(:, ind);
                vecR = datR(:, ind);
                
                [dGi, mGi] = getVecSum( vecG );
                mG = [mG, mGi/sum(abs(vecG))];

                [dRi, mRi] = getVecSum( vecR );
                mR = [mR, mRi/sum(abs(vecR))];
            end
            
            vR = abs(vR)
            vF = abs(vF)
            intG = intG(1:length(vR))
            intR = intR(1:length(vR))
            mG = mG(1:length(vR))
            mR = mR(1:length(vR))
            
            data{1} = [data{1} vR]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
            data{2} = [data{2} vF]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
            data{3} = [data{3} intG]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
            data{4} = [data{4} intR]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
            data{5} = [data{5} mG]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
            data{6} = [data{6} mR]%( mG > 0 & mR > 0 & mG < 3 & mR < 3)]
            
        end
    end
    
    %% Plot scatterplots
    name = strcat(dir, sprintf('fly%d_', i));
    for j=1:2
        if i == 1
            fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
        else
            fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
        end
        
        %iterate over green and red
        for k=1:2
            subplot(2, 2, 2*k-1)
            
            scatter(data{j}, data{2*k+1}, 3, 'filled')
            
            xlabel(names{j})
            ylabel(names{2*k+1})
            
            
            subplot(2, 2, 2*k)
            
            scatter(data{j}, data{2*k+2}, 3, 'filled')
            
            xlabel(names{j})
            ylabel(names{2*k+2})
            
        end         
        print(fig, strcat(name, strcat('scatter_', names{j})), '-dpdf');
        
        if i == 1
            box = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
        else
            box = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
        end
        %set(box, 'DefaultTextFontSize', 4);
        for k=1:4 %iterate over red and green
            subplot(4, 1, k)
            
            boxdat = [];
            labels = {};
            xlabs = [];
            distDat = [];%
            for m = 1:length(vels{j})-1
                if m == 1
                    dat = data{k+2}( data{j} == 0 );
                    dat0 = dat;
                    ps = [1];
                    groups = [ zeros(length(dat), 1) ];
                else
                    dat = data{k+2}( data{j} > vels{j}(m) & data{j} <= vels{j}(m+1) );
                    if length(dat) > 0
                        groups = vertcat( groups, (m-1)*ones(length(dat),1) );
                    end
                    
                end
                
                if length(dat) > 0
                    vertcat(distDat, dat);%
                    [r, p] = ttest2(dat, dat0, 'Tail', 'both', 'Vartype', 'unequal')
                    xlabs = [xlabs length(dat) p];
                    if k==1 || k==2
                        dat = dat / mean(dat0);
                    end
                    boxdat = [boxdat dat];
                    M = mean(dat);
                    dev = std(dat);
                    labels{m} = sprintf('%.2f m=%.2f(%.2f)', vels{j}(m+1), M, dev);
                end
            end

            xlab = ''
            for m = 1:length(ps)
                xlab = strcat(xlab, ',  n=%d p=%.2e') 
            end

            distributionPlot(transpose(boxdat), 'groups', transpose(groups), 'color', 'b',...
                'showMM', 5, 'xNames', labels)
            ylim( [0 max(boxdat)] )
            set(gca,'FontSize',8);
            %boxplot(boxdat, groups, 'Labels',labels, 'Whisker', 5)
            %set(findobj(gca,'Type','text'),'FontSize',18)
            xlabel( sprintf(xlab, xlabs), 'fontsize', 9)
            title(names{k+2}, 'FontSize', 16) 
            
        end
        box.PaperUnits = 'inches';
        box.PaperPosition = [0 0 8 11.5];
        print(box, strcat(name, strcat('violin_', names{j})), '-dpdf');
    end
end