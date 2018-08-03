function alignByVelPB(dir, green, red)
%given a directory and names of green and red channels, aligns one channel
%to the other (does it both ways) and plots the mean distributions with
%standard deviations. Uses alignPB. Puts lPB and rPB on same figure,
%parallel subfigures


try
    from_file = load(strcat(dir, 'cont'), 'alldata');
    alldata=from_file.alldata;
    
catch

    alldata = FlyDatLoad(2, 'PB');
    save(strcat(dir, 'cont'), 'alldata');
end

vels = [-20 -2*pi/3 -pi/3 -pi/6 0 0 pi/6 pi/3 2*pi/3 20]; %bin by rotational velocity

for i = 1:length(alldata{1}.allFlyData); %create new plot for each fly
    %% Get data
    fly = alldata{1}.allFlyData{i};
            
    %pool data for a given fly
    data = { {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]} };
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
   
        if length(trial) > 0 & max(trial.positionDatMatch.vF) > 0 %Only consider trials where the fly moves
            
            datG = flipud(trial.GROIaveMax)-1; %we're numbering from 1 at L9 to 18 at R9 on Tanya's diagram
            datR = flipud(trial.RROIaveMax)-1; %this means CCW rotation (+ve vrot) gives higher numbers
            %closed loop only
            vR = trial.positionDatMatch.vRot( trial.positionDatMatch.Closed(1:length(trial.positionDatMatch.vRot))== 1 );

            %should smoothen by SG
            smooth = 0;
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
                %consider each data point and add to appropriate binned
                %dataset (this is a pretty bad way of doing things)
                
                a = vels(k);
                b = vels(k+1);
                
                for m = 1:length(vR)
                
                    if a <= vR(m) & vR(m) <= b
                        dat = data(k);
                        data{k}{1} = [data{k}{1} datG(:,m)];
                        data{k}{2} = [data{k}{2} datR(:,m)];
                
                    end
                    
                end
                
            end
        end
    end
    
    
%% Plot with green aligned

    fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

    for m = 1:length(vels)-1
        sprintf('vel %.2f length %d', vels(m), length(data{m}{1}))
        
        [mGL, stdGL, mRL, stdRL, mGR, stdGR, mRR, stdRR, dirGL, magGL, dirRL, magRL, dirGR, magGR,...
            dirRR, magRR] = alignPB(data{m}{1}, data{m}{2}, green, red);                

        
        subplot(length(vels)-1, 2, 2*m-1)
        
        ymax = 1.05*max( [max(mGL(mGL~=0)+stdGL(mGL~=0)) max(mRL(mRL~=0)+stdRL(mRL~=0)) magGL magRL ] );
        ymin = 0.75*min( [min(mGL(mGL~=0)-stdGL(mGL~=0)) min(mRL(mRL~=0)-stdRL(mRL~=0)) magGL magRL ] );
                
        hold on
        x = [1:9,fliplr(1:9)];
        yy = [mGL-stdGL,fliplr(mGL+stdGL)];
        fill(x,yy,'g','facealpha',.1, 'LineStyle', 'none')
        
        yy = [mRL-stdRL,fliplr(mRL+stdRL)];
        fill(x,yy,'r','facealpha',.07, 'LineStyle', 'none')
        
        plot(1:9, mGL, 'g');
        plot(1:9, mRL, 'r');
        
        line([dirGL dirGL], [0 magGL], 'Color', 'g', 'LineWidth', 2)
        line([dirRL dirRL], [0 magRL], 'Color', 'r', 'LineWidth', 2)

        xlim([1 9])
        ylim([0 ymax])
        title( sprintf( strcat('L %.1f<vR<%.1f ', 'G ',...
            green, ' R ', red,' m%.2f std%.2f dG%.1f dR%.1f n%d'),...
            vels(m), vels(m+1), mGL(5), stdGL(5), dirGL, dirRL, length(data{m}{1}) ), 'FontSize', 8  )
        %legend({green, red})

%% Plot right
        subplot(length(vels)-1, 2, 2*m)
        
        ymax = 1.05*max( [max(mGR(mGR~=0)+stdGR(mGR~=0)) max(mRR(mRR~=0)+stdRR(mRR~=0)) magGR magRR ] );
        ymin = 0.75*min( [min(mGR(mGR~=0)-stdGR(mGR~=0)) min(mRR(mRR~=0)-stdRR(mRR~=0)) magGR magRR ] );
        hold on

        yy = [mGR-stdGR,fliplr(mGR+stdGR)];
        fill(x,yy,'g','facealpha',.1, 'LineStyle', 'none')
        
        yy = [mRR-stdRR,fliplr(mRR+stdRR)];
        fill(x,yy,'r','facealpha',.07, 'LineStyle', 'none')
        
        plot(1:9, mGR, 'g')
        plot(1:9, mRR, 'r')
        
        line([dirGR dirGR], [0 magGL], 'Color', 'g', 'LineWidth', 2)
        line([dirRR dirRR], [0 magRL], 'Color', 'r', 'LineWidth', 2)

        xlim([1 9])
        ylim([0 ymax])
        title( sprintf( strcat('R %.1f<vR<%.1f ', 'G ',...
            green, ' R ', red,' m%.2f std%.2f dG%.1f dR%.1f n%d'),...
            vels(m), vels(m+1), mGR(5), stdGR(5), dirGR, dirRR, length(data{m}{1}) ), 'FontSize', 8  )
        %legend({green, red})
    end
    
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 11.5];
    print(fig, strcat(dir, sprintf( 'fly%d_',i), green,'_aligned_', red ), '-dpdf');
    
%% Plot with red aligned    

    fig2 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
            
    for m = 1:length(vels)-1
        
        [mGL, stdGL, mRL, stdRL, mGR, stdGR, mRR, stdRR, dirGL, magGL, dirRL, magRL, dirGR, magGR,...
            dirRR, magRR] = alignPB(data{m}{2}, data{m}{1}, red, green);                
        
        subplot(length(vels)-1, 2, 2*m-1)
        
        ymax = 1.05*max( [max(mGL(mGL~=0)+stdGL(mGL~=0)) max(mRL(mRL~=0)+stdRL(mRL~=0)) magGL magRL ] );
        ymin = 0.75*min( [min(mGL(mGL~=0)-stdGL(mGL~=0)) min(mRL(mRL~=0)-stdRL(mRL~=0)) magGL magRL ] );
        
        hold on

        yy = [mGL-stdGL,fliplr(mGL+stdGL)];
        fill(x,yy,'g','facealpha',.1, 'LineStyle', 'none')
        
        yy = [mRL-stdRL,fliplr(mRL+stdRL)];
        fill(x,yy,'r','facealpha',.07, 'LineStyle', 'none')
        
        plot(1:9, mGL, 'g')
        plot(1:9, mRL, 'r')

        line([dirGL dirGL], [0 magGL], 'Color', 'g', 'LineWidth', 2)
        line([dirRL dirRL], [0 magRL], 'Color', 'r', 'LineWidth', 2)

        xlim([1 9])
        ylim([0 ymax])
        
        title( sprintf( strcat('L %.1f<vR<%.1f ', 'G ',...
            red, ' R ', green,' m%.2f std%.2f dG%.1f dR%.1f n%d'),...
            vels(m), vels(m+1), mGL(5), stdGL(5), dirGL, dirRL, length(data{m}{1}) ), 'FontSize', 8  )
        %legend({green, red})

%% Plot right        
        subplot(length(vels)-1, 2, 2*m)
        hold on
        
        ymax = 1.05*max( [max(mGR(mGR~=0)+stdGR(mGR~=0)) max(mRR(mRR~=0)+stdRR(mRR~=0)) magGR magRR ] );
        ymin = 0.75*min( [min(mGR(mGR~=0)-stdGR(mGR~=0)) min(mRR(mRR~=0)-stdRR(mRR~=0)) magGR magRR ] );

        yy = [mGR-stdGR,fliplr(mGR+stdGR)];
        fill(x,yy,'g','facealpha',.1, 'LineStyle', 'none')
        
        yy = [mRR-stdRR,fliplr(mRR+stdRR)];
        fill(x,yy,'r','facealpha',.07, 'LineStyle', 'none')
        
        plot(1:9, mGR, 'g')
        plot(1:9, mRR, 'r')
        
        line([dirGR dirGR], [0 magGR], 'Color', 'g', 'LineWidth', 2)
        line([dirRR dirRR], [0 magRR], 'Color', 'r', 'LineWidth', 2)

        xlim([1 9])
        ylim([0 ymax])
        title( sprintf( strcat('R %.1f<vR<%.1f ', 'G ',...
            red, ' R ', green,' m%.2f std%.2f dG%.1f dR%.1f n%d'),...
            vels(m), vels(m+1), mGR(5), stdGR(5), dirGR, dirRR, length(data{m}{1}) ), 'FontSize', 8  )
        %legend({green, red})
    end
    fig2.PaperUnits = 'inches';
    fig2.PaperPosition = [0 0 8 11.5];
    
    print(fig2, strcat(dir, sprintf( 'fly%d_',i), red,'_aligned_', green ), '-dpdf');
    
    
end




