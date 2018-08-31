function alignByVel1C(dir, cond)
%given a directory in which to save data and a container with data for an
%experiment, pools over flies and trials to align imaging data to max
%intensity (at n=5) and plots for vRot bins. Also quantifies FWHM and plots
%vs. vRot bins for cold and hot.

name = cond{1}.name;

vels = [-20 -2*pi/3 -pi/3 -pi/6 0 0 pi/6 pi/3 2*pi/3 20]; %pool according to velocity

pers = {'dark', 'CL', 'OL'};

conds = {'RT', '30C'};

for per = 1:3 %loop over periods of interest
    [datas, vRots] = align1C(cond, pers{per}); %get pooled aligned data

    fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    fwhms = {zeros(1, length(vels)-1), zeros(1, length(vels)-1)}; %initialize fwhm data structure
    for j = 1:2 %RT, 30C
        for m = 1:length(vels)-1 %iterate over bins
            a = vels(m);
            b = vels(m+1);
            
            select = (a <= vRots{j} & vRots{j} <= b); %find the indices corresponding to vRot of interest
            Gal = datas{j}(:,select); %get aligned intensity corresponding to bin
            
            mG = mean(Gal, 2)';
            stdG = zeros(1, 16);
            for i = 1:16
                stdG(i) = std(Gal(i,:));
            end

            [direcG, magG] = getVecSum(mG); %get PVA
            direcG = direcG * 16/(2*pi); %get direction in wedge units
            magG = magG / sum(abs(mG)); %get PVA magnitude
            
            [xs, HM] = FWHM(mG, 'min'); %find FWHM
            x1 = xs{1};
            x2 = xs{2};
            fwhm = x2-x1;
            
            fwhms{j}(m) = fwhm;

            ymax = 1.05*max( [max(mG+stdG) magG ] );
            ymin = 0.75*min( [min(mG-stdG) magG ] );

            subplot(length(vels)-1, 2, 2*(m-1)+j) %put RT, 30C next to each other

            hold on

            %plot standard deviation
            x = [1:16,fliplr(1:16)];
            yy = [mG-stdG,fliplr(mG+stdG)];
            fill(x,yy,'g','facealpha',.1, 'LineStyle', 'none')

            plot(1:16, mG, 'g')

            %plot mean
            line([direcG direcG], [0 magG], 'Color', 'g', 'LineWidth', 2)
            
            %plot FWHM
            line([x1 x1], [0 HM], 'Color', 'b', 'LineStyle', ':',  'LineWidth', 1)
            line([x2 x2], [0 HM], 'Color', 'b', 'LineStyle', ':', 'LineWidth', 1)
            line([x1 x2], [HM HM], 'Color', 'b', 'LineStyle', ':', 'LineWidth', 1)

            xlim([1 16])
            ylim([0 ymax])
            title( sprintf( strcat(name,conds{j},'-',pers{per},' %.1f:%.1f m=%.2f sd=%.2f d=%.2f n=%d FW=%.2f'),...
                vels(m), vels(m+1), mG(9), stdG(9), direcG, length(Gal), fwhm ), 'FontSize', 8  )
        end
    end

    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 11.5];
    print(fig, strcat( dir, 'aligned_', name, '_', pers{per} ), '-dpdf'); %save figure with all the plots
    
    fwhm_fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    %new figure for plotting how FWHM changes with vRot for RT and 30C
    xs = 1:length(vels)-1;
    
    subplot(2,1,1)
    plot(xs, fwhms{1}, 'bo--'); %RT is blue
    hold on
    plot(xs, fwhms{2}, 'ro--'); %30C is red
    xlim( [0, length(vels)] )

    
    %ylim( [ 0.8*min( min(fwhms{1}), min(fwhms{2}) ) 1.2*max( max(fwhms{1}), max(fwhms{2}) ) ] )
    
    ylim( [ 0.5*mean(fwhms{1}) 1.5*mean(fwhms{1}) ] ) %center the RT curve
    
    legend( {'RT', '30C'} )
    
    set(gca,'xtick',xs); 
    xlabs = {};
    for i = 1:length(vels)-1
        xlabs{i} = sprintf('%.1f:%.1f(%.0f)', vels(i), vels(i+1), 100*(fwhms{2}(i)-fwhms{1}(i))/fwhms{1}(i));
        %add percent change form cold to hot as label
    end
    set(gca,'xticklabel',xlabs, 'FontSize', 7);
    
    ylabel('FWHM', 'FontSize', 12);
    
    title( strcat(name,'-',pers{per}), 'FontSize', 15  )
    
    fwhm_fig.PaperUnits = 'inches';
    fwhm_fig.PaperPosition = [0 0 8 11.5];
    print(fwhm_fig, strcat( dir, 'fwhm_', name, '_', pers{per}, '_from_baseline' ), '-dpdf'); %save
    
end
