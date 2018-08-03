%go through directories with data and run analyses

clear

dirs = { '~/Documents/Imaging/Data_Dan/PB/PEN1_G_EPG_R_PB/',...
    '~/Documents/Imaging/Data_Dan/PB/PEN1_R_EPG_G_PB/',...
    '~/Documents/Imaging/Data_Dan/PB/PEN2_G_PEG_R_PB/',...
    '~/Documents/Imaging/Data_Dan/PB/PEN2_R_EPG_G_PB/',...
    '~/Documents/Imaging/Data_Dan/PB/D7_R_EPG_G_PB/',...
    '~/Documents/Imaging/Data_Dan/PB/D7_G_EPG_R_PB/' };

greens = { 'PEN1' 'EPG' 'PEN2' 'EPG' 'EPG' 'D7' };

reds = { 'EPG' 'PEN1' 'PEG' 'PEN2' 'D7' 'EPG' };

for i = 1:length(dirs)
    display(dirs{i})
    
    try
        from_file = load(strcat(dirs{i}, 'cont'), 'alldata');
        cond=from_file.alldata;    
    catch
        cond = FlyDatLoad(1);
        save(strcat(dirs{i}, 'cont'), 'alldata');
    end
    
    getDecayDelayEBPB(strcat(dirs{i}, 'decayDelay/'), cond, 7, 10, greens{i})
    %PB_analysis(dirs{i}, greens{i}, reds{i})
    %scatterVelPB(dirs{i}, greens{i}, reds{i})
    %alignByVelPB(dirs{i}, greens{i}, reds{i})
    
end
