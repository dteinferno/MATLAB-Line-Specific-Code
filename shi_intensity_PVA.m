%plots information on intensity (relative to RT), PVA and fwhm for each
%line in different velocity bins. Scatterplot with a point per trial. Also
%plots the variation of each quantity with velocity for a given line for
%both RT and 30C.

nWedge = 1

dirs = { '/Users/loaner/Documents/imaging/Data_Dan/shi/empty/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/EPG/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/GE/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEG/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEN2/'};

%EPG does not have much stationary info so crashes and we leave it out
dirs = { '/Users/loaner/Documents/imaging/Data_Dan/shi/empty/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEG/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEN2/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/GE/'};

vels = [0 0 pi/6 pi/3 2*pi/3 20]; %rotational velocities


for period = {'dark', 'CL', 'OL'}
    
    per = period{1}
    
    datas = cell(6, length(dirs)); %store numbers
    DFs = cell(6, length(dirs)); %store intensities
    names = cell(1,length(dirs));
    conds = cell(1,length(dirs));

    %% Get Data
    for d = 1:length(dirs)

        dir = dirs{d} %iterate over types

        try %load container
            from_file = load(strcat(dir, 'cond'), 'cond');
            cond=from_file.cond;    
        catch
            cond = FlyDatLoad(1);
            save(strcat(dir, 'cond'), 'cond');
        end
        
        conds{d} = cond;
        
        names{d} = cond{1}.name;

        for i = 1:length(cond{1}.allFlyData); %iterate over flies

            for k = 1:2 %iterate over trials. getTrialShi automatically adds 1 to trialID for 30C
                
                ind = 2*(i-1)+k;
                
                [~, DFs{ind, d}, datas{ind, d}] = getTrialShi(cond, per, i, k, 0, 0, nWedge); %get data and store for later
            end
        end
    end
                
    %% extract by vel
    

    intfig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    PVAfig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    FWHMfig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    
    figs = {intfig, PVAfig, FWHMfig};
    
    for m = 1:length(vels)-1
        
        ints = zeros(6, length(dirs)); %intensities
        PVAs = zeros(6, length(dirs)); %|PVA|
        FWHMs = zeros(6, length(dirs));
        
        a = vels(m);
        b = vels(m+1);
        
        for d = 1:length(dirs)        
            dir = dirs{d};
            for i = 1:length(conds{d}{1}.allFlyData); %iterate over flies
                for k = 1:2 %iterate over trials. getTrialShi automatically adds 1 to trialID for 30C
                    
                    ind = 2*(i-1)+k;
                    
                    data = datas{ind, d};
                    vRotsRT = abs(data{2});
                    vRots30 = abs(data{2+6});

                    DF = DFs{ind, d}{2}; %get 30C
                    size(DF)
                    
                    selectRT = (a <= vRotsRT & vRotsRT <= b); %indicies corresponding to bin of interest
                    select30 = (a <= vRots30 & vRots30 <= b);
                    
                    if any(selectRT ~=0) && any(select30 ~= 0) %if we have data
                        ints(ind, d) = mean(data{6+4}(select30)) / mean(data{4}(selectRT)); %normslize intensity to RT
                        PVAs(ind, d) = mean(data{6+5}(select30)); %dont normalize PVA
                        DF = DF(:,select30);
                        [means, stds, ~, ~] = alignEB(DF, DF);
                        [xs, HM] = FWHM(means, 'min'); %get half max from baseline
                        FWHMs(ind, d) = xs{2}-xs{1};
                        
                    end
                end
            end
        end
        
        
        %% Plot scatterplots
        
        plotdata = {ints, PVAs, FWHMs};
        
        for i = 1:3
            figure(figs{i})
            subplot(length(vels)-1, 2, m)
            hold on
            xlabs = '';
            val0 = plotdata{i}(:, 1);
            for d = 1:length(dirs)
                vals = plotdata{i}(:, d);
                vals
                vals = vals(vals~=0); %zeros only if fewer samples or no standing still; real data does not give zero exactly
                xs = linspace(-0.2, 0.2, length(vals)) + d;
                scatter(xs', vals, 'b');
                line([xs(1), xs(end)], [mean(vals), mean(vals)], 'Color', 'r', 'LineStyle', '-')
                line([xs(1), xs(end)], [mean(vals)-std(vals), mean(vals)-std(vals)], 'Color', 'r', 'LineStyle', ':')
                line([xs(1), xs(end)], [mean(vals)+std(vals), mean(vals)+std(vals)], 'Color', 'r', 'LineStyle', ':')
                [r, p] = ttest2(vals, val0(val0 ~= 0), 'Tail', 'both', 'Vartype', 'unequal');
                xlabs = strcat(xlabs, sprintf(' %.2f(%.2f):%.1e', mean(vals), std(vals), p));
            end  
            xlabel('samples')
            if i == 1
                ylabel('relative intensity')
                ylim( [0 1.5] )
            elseif i == 2
                ylabel('|PVA|')
                ylim( [0 1] )
            else
                ylabel('FWHM')
                ylim( [0 8] )
            end
            title( sprintf(strcat(' %.1f < |vRot| < %.1f ', per), a, b) )
            
            set(gca,'xtick',[1:length(dirs)]);
            set(gca,'xticklabel',names, 'FontSize', 9);
            
            xlabel( strcat(xlabs), 'fontsize', 7)
            
        end    
    end
    

    
    intfig.PaperUnits = 'inches';
    intfig.PaperPosition = [0 0 8 11.5];
    print(intfig, strcat( '/Users/loaner/Documents/imaging/Data_Dan/shi/plots/intensity_scatterplots_',...
        per,'_nWedge',num2str(nWedge)), '-dpdf');
    
    PVAfig.PaperUnits = 'inches';
    PVAfig.PaperPosition = [0 0 8 11.5];
    print(PVAfig, strcat( '/Users/loaner/Documents/imaging/Data_Dan/shi/plots/PVA_scatterplots_', per), '-dpdf');

    FWHMfig.PaperUnits = 'inches';
    FWHMfig.PaperPosition = [0 0 8 11.5];
    print(FWHMfig, strcat( '/Users/loaner/Documents/imaging/Data_Dan/shi/plots/FWHM_baseline_scatterplots_', per), '-dpdf');
    
    %% Plot change with velocity bin
    
    %%{
    for d = 1:length(dirs)
        [~, DFs, data, flyinds] = getDataShi(conds{d}, per, 0, 0, nWedge); %don't smoothen
        
        line_fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
        xs = 1:length(vels)-1;

        plotdata = { {zeros(1, length(vels)-1) zeros(1, length(vels)-1)} {zeros(1, length(vels)-1) zeros(1, length(vels)-1)} }; %contains int/PVA and cold/hot
        stddata = plotdata;
        
        for i = 1:2 %hot vs cold
            for m = 1:length(vels)-1
                a = vels(m);
                b = vels(m+1);
                vRots = abs(data{6*i-4});
                select = (a <= vRots & vRots <= b);
                plotdata{1}{i}(m) = mean(data{6*i-2}(select)); %intensity
                plotdata{2}{i}(m) = mean(data{6*i-1}(select)); %PVA magnitude
                stddata{1}{i}(m) = std(data{6*i-2}(select)); %intensity
                stddata{2}{i}(m) = std(data{6*i-1}(select)); %PVA magnitude                
            end
        end
        
        
        %% plot intensity
        subplot(2,1,1)
        hold on
        plot(xs, plotdata{1}{1}, 'bo--');
        plot(xs, plotdata{1}{2}, 'ro--');
        plot( [xs', xs'], [ (plotdata{1}{1}-stddata{1}{1})', (plotdata{1}{1}+stddata{1}{1})' ], 'b:' );
        plot( [xs', xs'], [ (plotdata{1}{2}-stddata{1}{2})', (plotdata{1}{2}+stddata{1}{2})' ], 'r:' );
        xlim( [0, length(vels)] );
        %ylim( [ 0.8*min( min(fwhms{1}), min(fwhms{2}) ) 1.2*max( max(fwhms{1}), max(fwhms{2}) ) ] )
        %ylim( [ 0.5*mean(plotdata{1}{1}) 1.5*mean(plotdata{1}{1}) ] )
        ylim( [ 0 2*mean(plotdata{1}{1}) ] );
        legend( {'RT', '30C'} );
        set(gca,'xtick',xs); 
        xlabs = {};
        for i = 1:length(vels)-1
            xlabs{i} = sprintf('%.1f:%.1f(%.0f)', vels(i), vels(i+1),...
                100*(plotdata{1}{2}(i)-plotdata{1}{1}(i))/plotdata{1}{1}(i)); %add percent change form cold to hot
        end
        set(gca,'xticklabel',xlabs, 'FontSize', 8);
        ylabel('Intensity', 'FontSize', 12);
        title( strcat(names{d},'-',per), 'FontSize', 15  );

        %% Plot PVA
        subplot(2,1,2)
        hold on
        plot(xs, plotdata{2}{1}, 'bo--');
        plot(xs, plotdata{2}{2}, 'ro--');
        plot( [xs', xs'], [ (plotdata{2}{1}-stddata{2}{1})', (plotdata{2}{1}+stddata{2}{1})' ], 'b:' );
        plot( [xs', xs'], [ (plotdata{2}{2}-stddata{2}{2})', (plotdata{2}{2}+stddata{2}{2})' ], 'r:' );
        xlim( [0, length(vels)] );
        ylim( [ 0 2*mean(plotdata{2}{1}) ] );
        legend( {'RT', '30C'} );
        set(gca,'xtick',xs); 
        xlabs = {};
        for i = 1:length(vels)-1
            xlabs{i} = sprintf('%.1f:%.1f(%.0f)', vels(i), vels(i+1),...
                100*(plotdata{2}{2}(i)-plotdata{2}{1}(i))/plotdata{2}{1}(i)); %add percent change form cold to hot
        end
        set(gca,'xticklabel',xlabs, 'FontSize', 8);
        ylabel('|PVA|', 'FontSize', 12);
        title( strcat(names{d},'-',per), 'FontSize', 15  ) 
        
        
        line_fig.PaperUnits = 'inches';
        line_fig.PaperPosition = [0 0 8 11.5];
        print(line_fig, strcat( dirs{d}, '/int_PVA_profiles/int_PVA_profile_', names{d},...
            '_', per, '_', 'nWedge', num2str(nWedge) ), '-dpdf');
        
        %error
    end
    %}
        
end
            
                    
