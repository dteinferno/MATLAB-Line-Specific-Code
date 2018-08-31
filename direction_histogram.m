%look at heading

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

conds = {'All' 'All_30C'};
j = 2

for periods = {'dark', 'CL', 'OL'}
    
    period = periods{1}
    
    %datas = cell(6, length(dirs)); %store numbers
    %DFs = cell(6, length(dirs)); %store intensities
    %names = cell(1,length(dirs));

    %% Get Data
    for d = 1:length(dirs)
        
        dirfig = figure('units','normalized','outerposition',[0 0 1 1])%, 'visible', 'off');

        dir = dirs{d} %iterate over types

        try %load container
            from_file = load(strcat(dir, 'cond'), 'cond');
            cond=from_file.cond;    
        catch
            cond = FlyDatLoad(1);
            save(strcat(dir, 'cond'), 'cond');
        end
        
        %conds{d} = cond;
        
        name = cond{1}.name;

        for i = 1:length(cond{1}.allFlyData); %iterate over flies

            for k = 2:3 %only consider last two warm trials 
                
                ind = 2*(i-1)+k-1;
                subplot(3,2,ind)
                smoothF = 5; %looking at directions, smooth a bit
                smoothV = 0; %don't smooth velocity
                
                [tPts,darkPer, OLPer, CLPer, CWPer, CCWPer, DF, heading, headingPlt, vRot, vF, stripePos, stripePosPlt, stripeJumps]...
                    = extractShiData(cond, i, conds{j}, k, smoothF, smoothV); %get data
                
                if strcmp(period, 'dark')
                    per = darkPer;
                elseif strcmp(period, 'OL')
                    per = OLPer;
                elseif strcmp(period, 'CL')
                    per = CLPer;
                else
                    display('per not recognized, please select dark or OL')
                end
                
                
                DF = DF(:, per);
                heading = heading(per);
                vF = vF(per);
                
                m = zeros(1, length(heading));
                direcs = zeros(1, length(heading));
                for ind = 1:length(heading);
                    vec = DF(:, ind);
                    [di, mi] = getVecSum( vec );
                    m(ind) = mi/sum(abs(vec)); %use normalized PVA amplitude
                    direcs(ind) = di-pi;
                end
                
                direcs = direcs( m > 0.25 & vF' ~= 0  );
                heading = heading( m > 0.25 & vF' ~= 0 );

                
                hold on
                h1 = histogram(direcs, 'BinEdges',linspace(-pi, pi, 17), 'Normalization', 'probability');
                h1.FaceColor = 'r'
                %h1(2).Color = 'r'
                h2 = histogram(heading, 'BinEdges',linspace(-pi, pi, 17), 'Normalization', 'probability');
                h2.FaceColor = 'b'
                %h2(2).Color = 'b'
                legend({'fluorescent directions', 'headings'})
                title( sprintf(strcat(name, '-', period, '-fly%d-', conds{j}, '-trial%d'), i, k) )
                xlim([-pi pi])
                ylim([0 0.4])
   
            end
        end
        
        dirfig.PaperUnits = 'inches';
        dirfig.PaperPosition = [0 0 8 11.5];
        print(dirfig, strcat( dir,'directions/histograms_', name, '_', period ), '-dpdf')
    end
end