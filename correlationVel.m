function rhos = correlationVel(dir, cond, offset)
%Given a directory to save files to and a container with data for an
%experiment, plots boxplots of rate of change of PVA vs vRot pooled across
%flies for each condition (OL, CL, dark). Also includes line of best fit
%and correlation. returns a list of R values.


name = cond{1}.name;
pers = { 'dark' 'OL' 'CL' };
cor = [];
conds = {'RT' '30C'};

rhos = cell(7,2,3);

for offset = [-3 -2 -1 0 1 2 3]

    for i = 1:2 %iterate over hot and cold. These get a figure each

        fig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
        
        
        for per = 1:3
            vRs = [];
            vDFs = [];
            for flyid = 1:cond{1}.numFlies
                for trialid = 1:2

                    [~, ~, data] = getTrialShi(cond, pers{per}, flyid, trialid, 9, 9, 16); %smoothen data since we're looking at directions

                    mag = data{6*i-1}; %PVA magnitudes

                    %tPts = data{6*i-5}(mag > 0.25); 
                    tPts = data{6*i-5};
                    dirs = data{6*i};
                    vR = data{6*i-4};
                    if offset > 0
                        vR = vR(1:end-offset);
                        dirs = dirs(offset+1:end);
                        mag = mag(offset+1:end);
                    elseif offset < 0
                        vR = vR(abs(offset)+1:end);
                        dirs = dirs(1:end-abs(offset));
                        mag = mag(1:end-abs(offset));
                    end

                    vR = vR(1:end-1);
                    mag = mag(1:end-1);

                    vDF = transpose( getRates(dirs, tPts, 2*pi) );

                    vRs = [vRs vR(mag > 0.25)]; %threshhold at PVA magnitude of 0.25 since we're looking at directions
                    vDFs = [vDFs vDF(mag > 0.25)];

                    rho = corr(transpose(vR(mag > 0.25)), transpose(vDF(mag > 0.25)))
                    
                    
                    rhos{offset+4, i, per} = [rhos{offset+4, i, per} rho];

                    %gets rate of change as difference between consecutive directions ovver difference in time


                    %vR = vR(abs(vDF) < 10); %consider everything else to be an outlier (this does not change results, only leaves out a few points)
                    %vDF = vDF(abs(vDF) < 10);

                    
                end
            end

            subplot(2,2,per)

            scatter(vRs, vDFs, 0.2, 'b')
            %plot rotation of fly vs bump
            hold on
            xlabel('vRot')
            ylabel('Rate of change of fluorescence direction')
            title( strcat( name, '-', pers{per},'-',conds{i} ) )

            p1 = polyfit(vRs,vDFs,1); %fit a line
            yfit1 = polyval(p1,vRs);
            rho = corr(transpose(vRs), transpose(vDFs)); %measure correlation

            %plot(vR, yfit1, 'k')
            plot([-10 10], [-10 10], 'k--')

            xlim([-10 10])
            ylim([-10,10])

            legend( {sprintf('R %.2f', rho)} )

            cor = [cor rho]; %return correlations in case they're needed

            
    

        end
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 8 11.5];
        print(fig, strcat(dir, 'scatter_corr_intensity_group_meeting_', name, '_', conds{i}, '_offset', num2str(offset)), '-dpdf'); %save
    end
end