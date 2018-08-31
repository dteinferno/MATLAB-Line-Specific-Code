
function jumpvel(dir, cond, smooth)

%Quantify rotational, foward velocity around a bar jump for all flies
%(Green fig2)
%affected by weight of ball?
%Get timeline' from -10 to 10s - for individual flies and averaged across
%flies
%consider maan turning velocity -2;0 and 0;2 for all flies, boxplot
%conside running a comparison between flies
%return list of pairwise changes in velocity


data = { [] [] [] [] [] [] [] [] [] [] [] [] };
DFs = {[] []};
conds = {'All' 'All_30C'};
names = {'RT', '30C'};
smoothF = 0;
smoothV = smooth;

name = cond{1}.name;


nsteps = 100;
npre = 15;

for j = 1:2; %iterate over conditions
    
    for jumpind = 1:2
    
        ts = []; %times for plotting
        vrots = [];
        vFs = [];
        pos = [];
        
        for i = 1:length(cond{1}.allFlyData); %iterate over flies

            for k = 1:2 %iterate over trials

                if j == 1
                    ind = k;
                else
                    ind = k+1; %don't use first trial from 30C
                end

                [tPts,darkPer, OLPer, CLPer, CWPer, CCWPer, DF, heading, headingPlt, vRot, vF, stripePos, stripePosPlt, stripeJumps]...
                    = extractShiData(cond, i, conds{j}, ind, smoothF, smoothV);

                jump = stripeJumps(jumpind);
                ts = [ts tPts(jump-nsteps+1:jump+nsteps)-tPts(jump)];
                vrots = [vrots vRot(jump-nsteps+1:jump+nsteps)];
                vFs = [vFs vF(jump-nsteps+1:jump+nsteps)];
                
                pos = [pos stripePos(jump-nsteps+1:jump+nsteps)-mean(stripePos(jump-nsteps-1:jump+nsteps))];
                
            end
        end

        newfig = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

        subplot(3,2,1)

        plot(ts, vrots)
        xlabel('time (s)')
        ylabel('vRot (rad/s)')
        title(strcat('Turning at jumps (trials) ', name, conds{j}))
        
        hold on
        plot([0, 0], [-pi/2, pi/2], 'k--')
        plot([ts(nsteps-npre,1), ts(nsteps-npre,1)], [-pi/2, pi/2], 'k:')
        plot([ts(nsteps+npre,1), ts(nsteps+npre,1)], [-pi/2, pi/2], 'k:')
        xlim( [ts(1) ts(end)] )

        subplot(3,2,2)
        plot(ts, vFs)
        xlabel('time (s)')
        ylabel('vF (units)')
        title(strcat('Forward movement at jumps (trials) ', name, names{j}))
        
        hold on
        plot([0, 0], [0, 1], 'k--')
        plot([ts(nsteps-npre,1), ts(nsteps-npre,1)], [0, 1], 'k:')
        plot([ts(nsteps+npre,1), ts(nsteps+npre,1)], [0, 1], 'k:')
        xlim( [ts(1) ts(end)] )
        
        
        subplot(3,2,3)
        plot(ts(:,1), mean(vrots, 2), 'b')
        xlabel('time (s)')
        ylabel('vRot (rad/s)')
        title(strcat('Mean turning at jumps ', name, names{j}))
        
        hold on
        %plot(ts(:,1), mean(pos, 2))
        plot([0, 0], [-pi/2, pi/2], 'r--')
        plot([ts(nsteps-npre,1), ts(nsteps-npre,1)], [-pi/2, pi/2], 'r:')
        plot([ts(nsteps+npre,1), ts(nsteps+npre,1)], [-pi/2, pi/2], 'r:')
        xlim( [ts(1) ts(end)] )
        
        
        subplot(3,2,4)
        plot(ts(:,1), mean(vFs, 2), 'b')
        hold on
        plot([0, 0], [0, 1], 'r--')
        plot([ts(nsteps-npre,1), ts(nsteps-npre,1)], [0, 1], 'r:')
        plot([ts(nsteps+npre,1), ts(nsteps+npre,1)], [0, 1], 'r:')
        
        xlabel('time (s)')
        ylabel('vF (units)')
        title(strcat('Mean forward movement at jumps ', name, names{j}))
        xlim( [ts(1), ts(end)] )

        subplot(3,2,5)

        pres = mean(vrots(nsteps-npre:nsteps, :), 1)
        posts = mean(vrots(nsteps+1:nsteps+1+npre, :), 1)
        
        [h, p] = ttest(pres, posts)

        xs = [linspace(0.9, 1.1, length(pres)) linspace(1.9, 2.1, length(posts))]

        scatter(xs, [pres posts])
        
        hold on
        plot([ones(1, length(pres)); 2*ones(1, length(pres))], [pres; posts], 'b:')
       
        ylabel('vRot (rad/s)')
        title(sprintf(strcat('dvRot (means(%.2f:%.2f) p=%.2e) ', name, names{j}), mean(pres), mean(posts), p ), 'FontSize', 9 )

        subplot(3,2,6)

        pres = mean(vFs(nsteps-npre:nsteps, :), 1);
        posts = mean(vFs(nsteps+1:nsteps+1+npre, :), 1);
        
        [h, p] = ttest(pres, posts)

        xs = [linspace(0.9, 1.1, length(pres)) linspace(1.9, 2.1, length(posts))];

        scatter(xs, [pres posts])
        hold on
        plot([ones(1, length(pres)); 2*ones(1, length(pres))], [pres; posts], 'b:')
        
        
        ylabel('vF (units)')
        title(sprintf(strcat('dvF (means(%.2f:%.2f) p=%.2e) ', name, names{j}), mean(pres), mean(posts), p ), 'FontSize', 9 )
        
        newfig.PaperUnits = 'inches';
        newfig.PaperPosition = [0 0 8 11.5];
        print(newfig, strcat(dir, strcat('jumpvel/jumpvel_', name, '_', names{j}, '_', num2str(jumpind))), '-dpdf');
    end
end