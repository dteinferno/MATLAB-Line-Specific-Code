function PB_analysis(dir, green, red)

cd ~/Documents/Imaging/Data_Dan

%dir = '~/Documents/Imaging/Data_Dan/PEN1_G_EPG_R_PB/';
%green = 'PEN1';
%red = 'EPG';

%dir = '~/Documents/Imaging/Data_Dan/PEN2_G_PEG_R_PB/';
%green = 'PEN2';
%red = 'PEG';

%dir = '~/Documents/Imaging/Data_Dan/PEN2_R_EPG_G_PB/';
%green = 'EPG';
%red = 'PEN2';

%dir = '~/Documents/Imaging/Data_Dan/D7_R_EPG_G_PB/';
%green = 'EPG';
%red = 'D7';

%dir = '~/Documents/Imaging/Data_Dan/D7_G_EPG_R_PB/';
%green = 'D7';
%red = 'EPG';

try
    from_file = load(strcat(dir, 'cont'), 'alldata');
    alldata=from_file.alldata;
    
catch

    alldata = FlyDatLoad(2, 'PB');
    save(strcat(dir, 'cont'), 'alldata');
end



ominusR = [] %offset green minus red at negative vrot (CW)
oplusR = []
ominusL = [] %offset green minus red at postitive crot (CCW)
oplusL = []
oAbsR = [] %the magnitude of the vector sum of green activity for |vrot| > pi/2 vs all vrot
oAbsL = []
actgR = []
actrR = []
actgL = []
actrL = []

slopes = []


for i = 1:length(alldata{1}.allFlyData);

    fly = alldata{1}.allFlyData{i};
    
    try
        L = length(fly.Dark);
    catch
        L = length(fly.All);
    end
    
    for j = L;

        try
            trial = fly.Dark{j};
        catch
            trial = fly.All{j};
        end
        
        if length(trial) > 0 & max(trial.positionDatMatch.vF) > 0 &...
                min(trial.positionDatMatch.vRot) < - 1.5*pi/2 &  max(trial.positionDatMatch.vRot) > 1.5*pi/2

            [a,b,c,d,e,f, g, h, k, l, m] = plot_PB_Kris(trial, green, red, i, j, dir);
            
            ominusR = [ominusR a] ;
            oplusR = [oplusR b];
            ominusL = [ominusL c]; 
            oplusL = [oplusL d];
            oAbsR = [oAbsR e];
            oAbsL = [oAbsL f];
            actgR = [actgR g];
            actrR = [actrR h];
            actgL = [actgL k];
            actrL = [actrL l];
            slopes = [slopes m];
            
        end
    end
end



        
boxes = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

subplot(2,2,1);

data = transpose([ominusR ; ominusL ; oplusR ; oplusL]);

m1 = mean(ominusR)
s1 = std(ominusR)
m2 = mean(ominusL)
s2 = std(ominusL)
m3 = mean(oplusR)
s3 = std(oplusR)
m4 = mean(oplusL)
s4 = std(oplusL)

boxplot(data, 'Labels',{ 'R-', 'L-', 'R+', 'L+'}, 'Whisker', 5)

xlabel(sprintf('%.2f+/-%.2f  %.2f+/-%.2f  %.2f+/-%.2f  %.2f+/-%.2f', m1,s1,m2,s2,m3,s3,m4,s4), 'fontsize', 8)

title('Offset')


subplot(2,2,2);

data = transpose([oAbsR ; oAbsL]);

m1 = mean(oAbsR)
s1 = std(oAbsR)
m2 = mean(oAbsL)
s2 = std(oAbsL)

boxplot(data, 'Labels',{sprintf('Right (%.2f+/-%.2f)',m1,s1),...
    sprintf('Left (%.2f+/-%.2f)',m2,s2)}, 'Whisker', 5)

title('|Offset|')


subplot(2,2,3);

data = transpose([actgR; actgL ; actrR ; actrL]);


m1 = mean(actgR)
s1 = std(actgR)
m2 = mean(actgL)
s2 = std(actgL)
m3 = mean(actrR)
s3 = std(actrR)
m4 = mean(actrL)
s4 = std(actrL)

boxplot(data, 'Labels',{ strcat(green, 'Right'), strcat(green, 'Left'),...
    strcat(red, 'Right'), strcat(red, 'Left')}, 'Whisker', 5)

xlabel(sprintf('%.2f+/-%.2f  %.2f+/-%.2f  %.2f+/-%.2f  %.2f+/-%.2f', m1,s1,m2,s2,m3,s3,m4,s4), 'fontsize', 8)

title('Activity')


subplot(2,2,4);

data = transpose([slopes( mod(1:length(slopes), 2)==0 ); slopes( mod(1:length(slopes), 2)==1 )]);

m1 = mean(data(:,1))
s1 = std(data(:,1))
m2 = mean(data(:,2))
s2 = std(data(:,2))

boxplot(data, 'Labels',{ sprintf('Right (%.2f+/-%.2f)',m1,s1),...
    sprintf('Left (%.2f+/-%.2f)',m2,s2)}, 'Whisker', 5)

title('Offset vs vRot slopes')


print(boxes, strcat(dir, 'boxplots_data'), '-dpdf');



save(strcat(dir, 'variables'),...
    'ominusR', 'oplusR', 'ominusL', 'oplusL', 'oAbsR', 'oAbsL', 'actgR', 'actrR', 'actgL', 'actrL', 'slopes')

