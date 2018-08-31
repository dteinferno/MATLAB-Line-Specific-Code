function rhos = correlationPos(dir, cond, offset)
%Given a directory to save files to and a container with data for an
%experiment, plots boxplots of direction of PVA vs fly heading pooled across
%flies for each condition (OL, CL, dark). Also includes line of best fit
%and correlation. returns a list of R values.


name = cond{1}.name;
pers = { 'dark' 'OL' 'CL' };
cor = [];
conds = {'RT' '30C'};


rhos = cell(7,2,3);

for offset = [-3 -2 -1 0 1 2 3]

    for i = 1:2 %iterate over hot and cold. These get a figure each


        for flyid = 1:3
            for trialid = 1:2
            fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
                for per = 1:3

                    [headings, ~, data] = getTrialShi(cond, pers{per}, flyid, trialid, 7, 7, 16); %smoothen data since we're looking at directions

                    mag = data{6*i-1}; %PVA magnitudes

                    %tPts = data{6*i-5}(mag > 0.25); 
                    tPts = data{6*i-5};
                    dirs = data{6*i}-pi;
                    headings = headings{i}';

                    if offset > 0
                        headings = headings(1:end-offset);
                        dirs = dirs(offset+1:end);
                        mag = mag(offset+1:end);
                    elseif offset < 0
                        headings = headings(abs(offset)+1:end);
                        dirs = dirs(1:end-abs(offset));
                        mag = mag(1:end-abs(offset));
                    end

                    difs = getOffset(dirs, headings, 2*pi); %get closest distance to posts from pres
                    dirs = headings+difs; %this allows us to 'wrap around' for the regression

                    %Headings = [Headings headings(mag > 0.25)]; %threshhold at PVA magnitude of 0.25 since we're looking at directions
                    %Poss = [Poss dirs(mag > 0.25)];

                    subplot(2,2,per)

                    scatter(headings, dirs, 'b')
                    %plot rotation of fly vs bump
                    hold on
                    xlabel('Heading')
                    ylabel('PVA direction')
                    title( strcat( name, '-', pers{per},'-',conds{i} ) )

                    p1 = polyfit(headings,dirs,1); %fit a line
                    yfit1 = polyval(p1,headings);
                    rho = corr(transpose(headings), transpose(dirs)); %measure correlation

                    %plot(vR, yfit1, 'k')
                    plot([-pi pi], [-pi pi], 'k--')

                    xlim([-pi pi])
                    ylim([-2*pi,2*pi])

                    legend( {sprintf('Slope %.2f R %.2f', p1(1), rho)} )

                    cor = [cor rho]; %return correlations in case they're needed

                    rhos{offset+4, i, per} = [rhos{offset+4, i, per} rho];

                    
                end
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0 0 8 11.5];
                print(fig, strcat(dir, 'scatter_corr_intensity_', name, '_', conds{i},...
                    num2str(flyid), ',', num2str(trialid), '_offset', num2str(offset)), '-dpdf'); %save
            end

        end
        
    end
end