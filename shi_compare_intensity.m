
%Plots a comparison diagram for GCaMP intensities and |PVA| between
%different shibire samples at 30C
clear

nWedge = 1

dirs = { '/Users/loaner/Documents/imaging/Data_Dan/shi/empty/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/EPG/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/GE/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEG/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEN2/'};


per = 'dark'; %specify per manually
datas = {[] [] [] []};

ind = 0;
for i = [1 4 5 3] %iterate over directories of interest
    ind = ind + 1; %have an index that counts non-empty data
    dir = dirs{i}
    
    try%get data if already stored
        from_file = load(strcat(dir, 'cond'), 'cond');
        cond=from_file.cond;    
    catch
        cond = FlyDatLoad(1);
        save(strcat(dir, 'cond'), 'cond');
    end
    
    data = scatterVelEB1C(strcat(dir, 'intensities/'), cond, per, nWedge); %get data by running scatterVelEB1C; also plots figures itself...
    
    datas{ind} = data; %4,5,6 contains our vRot, int, PCA. store by index
end

fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');


for j = 1:2 %run over intensity, PCA. These get a figure each. Normalize intensity
    subplot(3, 1, j)

    if j == 1
        ymax = 5.0
        %ymax = 6.0
    else
        ymax = 1.0
    end
    
    boxdat = [];
    labels = {};
    xlabs = [];
    distDat = [];%
    colmap = {};
    
    for k=1:4 %run over empty, PEG, PEN2, GE

        dat1 = datas{k}{4+j}( datas{k}{4} < pi/1000 ); 
        dat2 = datas{k}{4+j}( datas{k}{4} > pi/5 ); %compare standing completely still to  moving faster than pi/5
        
        if k == 1
            dat01 = dat1
            dat02 = dat2
            groups =  vertcat(zeros(length(dat1), 1), ones(length(dat2),1));
        else
            groups = vertcat( groups, (2*k-2)*ones(length(dat1),1) );
            groups = vertcat( groups, (2*k-1)*ones(length(dat2),1) );
        end

        if length(dat1) > 0
            vertcat(distDat, dat1);%
            vertcat(distDat, dat2);
                       
            if j == 1 %intensities
                [r1, p1] = ttest2(dat1, dat01, 'Tail', 'left', 'Vartype', 'unequal') %check IF LOWER than control
                [r2, p2] = ttest2(dat2/mean(dat1), dat02/mean(dat01), 'Tail', 'both', 'Vartype', 'unequal') %compare normalized
            else
                [r1, p1] = ttest2(dat1, dat01, 'Tail', 'both', 'Vartype', 'unequal') %check if different from control
                [r2, p2] = ttest2(dat2, dat02, 'Tail', 'both', 'Vartype', 'unequal')
            end
            xlabs = [xlabs length(dat1) p1 length(dat2) p2];

            if j == 1
                %dat1 = dat1 / mean(dat01);
                %dat2 = dat2 / mean(dat01);
                
                dat2 = dat2 / mean(dat1);
                dat1 = dat1 / mean(dat1);    
            end

            
            boxdat = [boxdat dat1 dat2];
            M1 = mean(dat1);
            dev1 = std(dat1);
            M2 = mean(dat2);
            dev2 = std(dat2);
            labels{2*k-1} = sprintf('%.2f(%.2f)', M1, dev1);
            labels{2*k} = sprintf('%.2f(%.2f)', M2, dev2);
            colmap{2*k-1} = 'b';
            colmap{2*k} = [0 1 1];
        end
    end

    xlab = '';
    for m = 1:8
        xlab = strcat(xlab, ',  n=%d p=%.2e');
    end

    distributionPlot(transpose(boxdat), 'groups', transpose(groups), 'color', colmap,...
        'showMM', 5, 'xNames', labels)
    ylim( [0 ymax] )
    set(gca,'FontSize',6);
    %boxplot(boxdat, groups, 'Labels',labels, 'Whisker', 5)
    %set(findobj(gca,'Type','text'),'FontSize',18)
    xlabel( sprintf(xlab, xlabs), 'fontsize', 6)
    if j == 1
        title('Intensity empty PEG PEN2 GE (vRot = 0.0, > pi/5)', 'FontSize', 16)
    else
        title('|PVA| empty PEG PEN2 GE (vRot = 0.0, > pi/5)', 'FontSize', 16)
    end
end
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 11.5];
print(fig, strcat('/Users/loaner/Documents/imaging/Data_Dan/shi/plots/compare_intensity_PEG_PEN2_GE', per, '_nWedge', num2str(nWedge)), '-dpdf');