dirs = { '/Users/loaner/Documents/imaging/Data_Dan/shi/empty/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/GE/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEG/',...
    '/Users/loaner/Documents/imaging/Data_Dan/shi/PEN2/'};

%{

pers = { 'dark' 'OL' 'CL' };
cor = [];
conds = {'RT' '30C'};
names = {'Ctrl', 'GE', 'PEG', 'PEN2'};

plotdir = '/Users/loaner/Documents/imaging/Data_Dan/shi/plots/';

rhos = cell(1, length(dirs));
for i = 1:length(dirs);
    dir = dirs{i}
    
    try
        from_file = load(strcat(dir, 'cond'), 'cond');
        cond=from_file.cond;    
    catch
        cond = FlyDatLoad(1);
        save(strcat(dir, 'cond'), 'cond');
    end
    
    rhos{i} = correlationVel(strcat(dir, 'correlation_vel/'), cond, 0);
end
    
%}



for i = 1:7 %offsets
    for j = 1:2 %hot vs cold

        scatterfig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
        for per = 1:3 %period

            subplot(3,1,per)
            hold on
            xlabs = '';
            vals0 = rhos{1}{i,j,per};

            for trial = 1:length(dirs)
                display(strcat('offset ', num2str(i-4), ' condition', num2str(j), ' period', num2str(per), ' trial', num2str(trial)))
                vals = rhos{trial}{i,j,per}
                xs = linspace(-0.2, 0.2, length(vals)) + trial;
                scatter(xs, vals, 'b');
                line([xs(1), xs(end)], [mean(vals), mean(vals)], 'Color', 'r', 'LineStyle', '-')
                line([xs(1), xs(end)], [mean(vals)-std(vals), mean(vals)-std(vals)], 'Color', 'r', 'LineStyle', ':')
                line([xs(1), xs(end)], [mean(vals)+std(vals), mean(vals)+std(vals)], 'Color', 'r', 'LineStyle', ':')
                [r, p] = ttest2(vals, vals0, 'Tail', 'both', 'Vartype', 'unequal');
                xlabs = strcat(xlabs, sprintf(' %.2f(%.2f):%.1e', mean(vals), std(vals), p));

            end
            title( pers{per} )
            xlim( [0 length(dirs)+1] )


            set(gca,'xtick',1:length(dirs));
            set(gca,'xticklabel',names, 'FontSize', 9);

            xlabel( strcat(xlabs), 'fontsize', 7)
            ylabel('vRot vs. PVA rate of change correlation')


        end

        scatterfig.PaperUnits = 'inches';
        scatterfig.PaperPosition = [0 0 8 11.5];
        print(scatterfig, strcat( plotdir, 'scattercor_vel_offset',num2str(i-4), '_', conds{j} ), '-dpdf');


    end
end