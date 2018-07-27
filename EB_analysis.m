function EB_analysis(dir, green, red)

%dir = '~/Documents/Imaging/Data_Dan/PEN2_G_PEG_R_EB/';
%green = 'PEN2';
%red = 'PEG';

%dir = '~/Documents/Imaging/Data_Dan/PEN1_G_EPG_R_EB/';
%green = 'PEN1';
%red = 'EPG';

%dir = '~/Documents/Imaging/Data_Dan/EPG_G_GE_R_EB/';
%green = 'EPG';
%red = 'GE';

%dir = '~/Documents/Imaging/Data_Dan/PEN2_G_EPG_R_EB/';
%green = 'PEN2';
%red = 'EPG';

%dir = '~/Documents/Imaging/Data_Dan/PEN2_R_EPG_G_EB/';
%green = 'EPG';
%red = 'PEN2';

try
    from_file = load(strcat(dir, 'cont'), 'alldata');
    alldata=from_file.alldata;
    
catch

    alldata = FlyDatLoad(2, 'EB');
    save(strcat(dir, 'cont'), 'alldata');
end

offset_neg = [] %offset green minus red at negative vrot (CW)
std_neg = []
offset_pos = [] %offset green minus red at postitive crot (CCW)
std_pos = []
actG = [] %the magnitude of the vector sum of green activity for |vrot| > pi/2 vs all vrot
actR = []
slopes = []


for i = 1:length(alldata{1}.allFlyData);

    fly = alldata{1}.allFlyData{i};
    
    try
        L = length(fly.Dark)
    catch
        L = length(fly.All)
    end
    
    for j = 1:L;
        
        try
            trial = fly.Dark{j}
        catch
            trial = fly.All{j}
        end
        
        if length(trial) > 0 & max(trial.positionDatMatch.vF) > 0

            [a,b,c,d,e,f g] = plot_vel_Kris(trial, green, red, i, j, dir)

            offset_neg = [offset_neg a];
            std_neg = [std_neg b];
            offset_pos = [offset_pos c];
            std_pos = [std_pos d];
            actG = [actG e];
            actR = [actR f];
            slopes = [slopes g];
        end
    end
end

        
boxes = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');

subplot(2,3,1);

data = transpose([offset_neg ; offset_pos]);

boxplot(data, 'Labels',{'Negative','Positive'}, 'Whisker', 5)

title('Offset')


subplot(2,3,2);

data = transpose([std_neg ; std_pos]);

boxplot(data, 'Labels',{'Negative','Positive'}, 'Whisker', 5)

title('Std')

subplot(2,3,3);

data = transpose([actG ; actR]);

boxplot(data, 'Labels',{green,red}, 'Whisker', 5)
mean(actG)
std(actG)
mean(actR)
std(actR)

title('Activity')

print(boxes, strcat(dir, 'boxplots_data'), '-dpdf');

save(strcat(dir, 'variables'),...
    'offset_neg', 'offset_pos', 'std_neg' ,'std_pos', 'actG', 'actR', 'slopes')




