function cor = correlationVel(dir, cond)
%Given a directory to save files to and a container with data for an
%experiment, plots boxplots of rate of change of PVA vs vRot pooled across
%flies for each condition (OL, CL, dark). Also includes line of best fit
%and correlation. returns a list of R values.


name = cond{1}.name;
pers = { 'dark' 'OL' 'CL' };
cor = [];
conds = {'RT' '30C'};

for i = 1:2 %iterate over hot and cold. These get a figure each

    fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    for per = 1:3 %iterate over periods

        data = getDataShi(cond, pers{per}, 7, 7); %smoothen data since we're looking at directions

        mag = data{6*i-1}; %PVA magnitudes

        tPts = data{6*i-5}(mag > 0.25); %threshhold at PVA amgnitude of 0.25 since we're looking at directions
        vR = data{6*i-4}(mag > 0.25);
        dirs = data{6*i}(mag > 0.25);

        vDF = transpose( getRates(dirs, tPts, 2*pi) );
        %gets rate of change as difference between consecutive directions ovver difference in time

        vR = vR(abs(vDF) < 10); %consider everything else to be an outlier (this does not change results, only leaves out a few points)
        vDF = vDF(abs(vDF) < 10);

        subplot(2,2,per)

        scatter(vR, vDF, 'b')
        %plot rotation of fly vs bump
        hold on
        xlabel('vRot')
        ylabel('Rate of change of fluorescence direction')
        title( strcat( name, '-', pers{per},'-',conds{i} ) )

        p1 = polyfit(vR,vDF,1); %fit a line
        yfit1 = polyval(p1,vR);
        rho = corr(transpose(vR), transpose(vDF)); %measure correlation

        plot(vR, yfit1, 'k')
        legend( {sprintf('Slope %.2f R %.2f', p1(1), rho)} )

        cor = [cor rho]; %return correlations in case they're needed
    end

    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 11.5];
    print(fig, strcat(dir, 'scatter_corr_intensity_', name, '_', conds{i}), '-dpdf'); %save
end