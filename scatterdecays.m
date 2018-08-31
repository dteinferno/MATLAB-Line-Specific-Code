%function scatterdecays(plotdir, region, nWedge)

nWedge = 3

%region = 'PB'
%plotdir = '~/Documents/Imaging/Data_Dan/PB/'

%region = 'shi';
%plotdir = '~/Documents/Imaging/Data_Dan/shi/plots/';

region = 'EB'
plotdir = '~/Documents/Imaging/Data_Dan/EB/'

pers = { 'dark', 'OL', 'CL' };
taus = { {} {} {} };
conds = {'RT', '30C'}

if strcmp(region, 'PB')
    dirs = { '~/Documents/Imaging/Data_Dan/PB/PEN2_R_EPG_G_PB/',...
        '~/Documents/Imaging/Data_Dan/PB/PEN1_R_EPG_G_PB/',...
    '~/Documents/Imaging/Data_Dan/PB/D7_R_EPG_G_PB/',...
    '~/Documents/Imaging/Data_Dan/PB/PEN1_G_EPG_R_PB/',...
    '~/Documents/Imaging/Data_Dan/PB/PEN2_G_PEG_R_PB/',...
    '~/Documents/Imaging/Data_Dan/1C/PEG_PB/',...
    '~/Documents/Imaging/Data_Dan/PB/D7_G_EPG_R_PB/'};

    names = { 'EPG' 'EPG' 'EPG' 'PEN1' 'PEN2' 'PEG' 'D7'};
    
elseif strcmp(region, 'EB')
    dirs = {'~/Documents/Imaging/Data_Dan/EB/PEN1_R_EPG_G_EB/',...
        '~/Documents/Imaging/Data_Dan/EB/EPG_G_GE_R_EB/',...
        '~/Documents/Imaging/Data_Dan/EB/PEN2_R_EPG_G_EB/',...
        '~/Documents/Imaging/Data_Dan/EB/PEN2_G_EPG_R_EB/',...
        '~/Documents/Imaging/Data_Dan/EB/PEN2_G_PEG_R_EB/',...
        '~/Documents/Imaging/Data_Dan/EB/PEN1_G_EPG_R_EB/',...
        '~/Documents/Imaging/Data_Dan/1C/PEG_EB/',...
        '~/Documents/Imaging/Data_Dan/EB/EPG_R_GE_G_EB/'};

    names = {'EPG' 'EPG' 'EPG' 'PEN2' 'PEN2' 'PEN1' 'PEG' 'GE'};
    
elseif strcmp(region, 'shi')
    dirs = { '/Users/loaner/Documents/imaging/Data_Dan/shi/empty/',...
        '/Users/loaner/Documents/imaging/Data_Dan/shi/EPG/',...
        '/Users/loaner/Documents/imaging/Data_Dan/shi/GE/',...
        '/Users/loaner/Documents/imaging/Data_Dan/shi/PEG/',...
        '/Users/loaner/Documents/imaging/Data_Dan/shi/PEN2/'};
    
    names = {'empty' 'EPG' 'GE' 'PEG' 'PEN2'};
    taus = { {{} {} {}}  {{} {} {}} };
else
    display('region not recognized')
end



for i = 1:length(dirs);
    dir = dirs{i}
    %dir = '~/Documents/Imaging/Data_Dan/1C/PEG_PB/'
    %dir = '~/Documents/Imaging/Data_Dan/1C/PEG_EB/'
    
    if strcmp(region, 'shi') || strcmp('PEG_', dir(end-6:end-3))
        try
            from_file = load(strcat(dir, 'cond'), 'cond');
            cond=from_file.cond;    
        catch
            cond = FlyDatLoad(1);
            save(strcat(dir, 'cond'), 'cond');
        end
    else
        try
            from_file = load(strcat(dir, 'cont'), 'alldata');
            cond=from_file.alldata;    
        catch
            cond = FlyDatLoad(1);
            save(strcat(dir, 'cont'), 'alldata');
        end
    end
    
    if strcmp(region, 'PB')
        [newtaus, tPts] = getDecayDelayPB(strcat(dirs{i}, 'decayDelay/'), cond, 7, 10, names{i}, nWedge);
        for per = pers
            nts = newtaus(per{1});
            if ~isempty(nts)
                newtaus(per{1}) = [ nts(:,1) ; nts(:,2)  ]; %combine left and right data
            end
        end
            
    elseif strcmp(region, 'EB')
        newtaus = getDecayDelayEB(strcat(dirs{i}, 'decayDelay/'), cond, 7, 10, names{i}, nWedge);

    elseif strcmp(region, 'shi')
        newtaus = getDecayDelay(strcat(dirs{i}, 'decayDelay/'), cond, 3, 10, nWedge);
        
    else
        display('region not recognized')
    end
    
    for per = 1:length(pers)
        if strcmp(region, 'shi')
            dat1 = newtaus{1}(pers{per});    
            dat1 = dat1( dat1 < mean(dat1)+2.5*std(dat1) & dat1 > mean(dat1)-2.5*std(dat1) ); %remove outliers to avoid distorting data
            taus{1}{per}{i} = dat1;
            dat2 = newtaus{2}(pers{per});    
            dat2 = dat2( dat2 < mean(dat2)+2.5*std(dat2) & dat2 > mean(dat2)-2.5*std(dat2) ); %remove outliers to avoid distorting data
            taus{2}{per}{i} = dat2;
        else
            dat = newtaus(pers{per});    
            dat = dat( dat < mean(dat)+2.5*std(dat) & dat > mean(dat)-2.5*std(dat) ); %remove outliers to avoid distorting data
            taus{per}{i} = dat;
        end
    end
end

if ~strcmp(region, 'shi')
    taus = {taus}
end

for t = 1:length(taus)
    scatterfig = figure('units','normalized','outerposition',[0 0 1 1])%, 'visible', 'off');
    for per = 1:3
        subplot(3, 1, per)
        hold on
        xlabs = '';
        val0 = [];
        for i = 1:length(dirs)
            if isempty(val0)
                val0 = taus{t}{per}{i}; %first group with entries
            end

            vals = taus{t}{per}{i};
            if ~isempty(vals)
                xs = linspace(-0.2, 0.2, length(vals)) + i;
                scatter(xs', vals, 'b');
                line([xs(1), xs(end)], [mean(vals), mean(vals)], 'Color', 'r', 'LineStyle', '-')
                line([xs(1), xs(end)], [mean(vals)-std(vals), mean(vals)-std(vals)], 'Color', 'r', 'LineStyle', ':')
                line([xs(1), xs(end)], [mean(vals)+std(vals), mean(vals)+std(vals)], 'Color', 'r', 'LineStyle', ':')
                [r, p] = ttest2(vals, val0, 'Tail', 'both', 'Vartype', 'unequal');
                xlabs = strcat(xlabs, sprintf(' %.2f(%.2f):%.1e', mean(vals), std(vals), p));
            end
        end

        title( pers{per} )
        xlim( [0 length(dirs)+1] )
        
        %{
        if strcmp(region, 'PB')
            ylim( [-10 20] )
        elseif strcmp(region, 'EB')
            ylim( [0 10] )
        else
            ylim( [-5 15] )
        end
        %}  

        set(gca,'xtick',1:length(dirs));
        set(gca,'xticklabel',names, 'FontSize', 9);

        xlabel( strcat(xlabs), 'fontsize', 7)
        ylabel('Rate constants for decay (s^(-1))')
    end

    scatterfig.PaperUnits = 'inches';
    scatterfig.PaperPosition = [0 0 8 11.5];
    print(scatterfig, strcat( plotdir, 'scatterdecay_',region, '_', conds{t}, '_','nWedge',num2str(nWedge) ), '-dpdf');
end
     
        
      