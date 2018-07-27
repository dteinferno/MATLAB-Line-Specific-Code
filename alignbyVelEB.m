function alignbyVelEB(dir, green, red)

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



vels = [-20 -2*pi/3 -pi/3 -pi/6 0 0 pi/6 pi/3 2*pi/3 20];

for i = 1:length(alldata{1}.allFlyData);

    fly = alldata{1}.allFlyData{i};
            
    data = { {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]} }
    
    try
        L = length(fly.Dark)
    catch
        L = length(fly.All)
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


            smooth = 3
            if smooth > 0
                s = size(datG);
                s = s(1);
                display('smoothening data');
                vR = Smooth(vR, smooth);
                newG = [];
                newR = [];
                for ind = 1:s;
                    g = datG(ind,:);
                    newG(ind,:) = Smooth(g, smooth);
                    r = datR(ind,:);
                    newR(ind,:) = Smooth(r, smooth);
                end
                datG = newG;
                datR = newR;    
            end            
            

            for k = 1:length(vels)-1
                
                a = vels(k);
                b = vels(k+1);
                
                for m = 1:length(vR)
                
                    if a <= vR(m) & vR(m) <= b
                        dat = data(k)
                        data{k}{1} = [data{k}{1} datG(:,m)];
                        data{k}{2} = [data{k}{2} datR(:,m)];
                
                    end
                    
                end
                
            end
        end
    end
    
    if i == 1
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
    else
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
    end
    
    for m = 1:length(vels)-1
        sprintf('vel %.2f length %d', vels(m), length(data{m}{1}))
        
        [mG, stdG, mR, stdR] = alignEB(data{m}{1}, data{m}{2}, green, red);

        
        [direcG, magG] = getVecSum(mG); 
        [direcR, magR] = getVecSum(mR);
        direcG = direcG * 16/(2*pi);
        direcR = direcR * 16/(2*pi);
        magR = magR / sum(mR) * 2; %get 15 times PVA
        magG = magG / sum(mG) * 2;
                
        ymax = 1.05*max( [max(mG+stdG) max(mR+stdR) magG magR ] );
        ymin = 0.75*min( [min(mG-stdG) min(mR-stdR) magG magR ] );
        
        subplot(length(vels)-1, 1, m)
                
        hold on

        x = [1:16,fliplr(1:16)];
        yy = [mG-stdG,fliplr(mG+stdG)];
        fill(x,yy,'g','facealpha',.1, 'LineStyle', 'none')
        
        yy = [mR-stdR,fliplr(mR+stdR)];
        fill(x,yy,'r','facealpha',.07, 'LineStyle', 'none')
        
        plot(1:16, mG, 'g')
        plot(1:16, mR, 'r')
        
        line([direcG direcG], [0 magG], 'Color', 'g', 'LineWidth', 2)
        line([direcR direcR], [0 magR], 'Color', 'r', 'LineWidth', 2)

        xlim([1 16])
        ylim([0 ymax])
        title( sprintf( strcat('EB %.2f<vR<%.2f ', 'G ',...
            green, ' R ', red,' m%.2f std%.2f dG%.2f dR%.2f n%d'),...
            vels(m), vels(m+1), mG(9), stdG(9), direcG, direcR, length(data{m}{1}) ), 'FontSize', 8  )
        %legend({green, red})

    end
    
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 11.5];
    print(fig, strcat(dir, sprintf( 'fly%d_',i), green,'_aligned_', red ), '-dpdf');

%% Align to red
    if i == 1
        fig2 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
    else
        fig2 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
    end
        
    for m = 1:length(vels)-1
        
        [mG, stdG, mR, stdR ] = alignEB(data{m}{2}, data{m}{1}, green, red); 
        
        [direcG, magG] = getVecSum(mG); 
        [direcR, magR] = getVecSum(mR);
        direcG = direcG * 16/(2*pi);
        direcR = direcR * 16/(2*pi);
        magR = magR / sum(mR) * 2;
        magG = magG / sum(mG) * 2;
        
        ymax = 1.05*max( [max(mG+stdG) max(mR+stdR) magG magR ] );
        ymin = 0.75*min( [min(mG-stdG) min(mR-stdR) magG magR ] );
        
        subplot(length(vels)-1, 1, m)
                
        hold on

        yy = [mG-stdG,fliplr(mG+stdG)];
        fill(x,yy,'g','facealpha',.1, 'LineStyle', 'none')
        
        yy = [mR-stdR,fliplr(mR+stdR)];
        fill(x,yy,'r','facealpha',.07, 'LineStyle', 'none')
        
        plot(1:16, mG, 'g')
        plot(1:16, mR, 'r')
        
        line([direcG direcG], [0 magG], 'Color', 'g', 'LineWidth', 2)
        line([direcR direcR], [0 magR], 'Color', 'r', 'LineWidth', 2)

        xlim([1 16])
        ylim([0 ymax])
        title( sprintf( strcat('EB %.2f<vR<%.2f ', 'G ',...
            red, ' R ', green,' m%.2f std%.2f dG%.2f dR%.2f n%d'),...
            vels(m), vels(m+1), mG(9), stdG(9), direcG, direcR, length(data{m}{1}) ), 'FontSize', 8   )
        %legend({green, red})

    end
    fig2.PaperUnits = 'inches';
    fig2.PaperPosition = [0 0 8 11.5];
    
    print(fig2, strcat(dir, sprintf( 'fly%d_',i), red,'_aligned_', green ), '-dpdf');
    
    
end